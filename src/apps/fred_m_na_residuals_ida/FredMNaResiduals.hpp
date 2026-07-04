#pragma once
#include "FredMNaState.hpp"
#include "FredMNaStressStrain.hpp"
#include "fuelpelletmaterial/UPuZr.hpp"
#include "claddingmaterial/HT9.hpp"
#include "platform/RodResiduals.hpp"

namespace fred {

constexpr double MNA_BUP_SCALE = 8.64e10; // J per MWd/kgU

// -----------------------------------------------------------------------
// FredMNaResiduals — DAE residual assembler for FRED-M-Na.
//
// Inherits shared infrastructure from RodResiduals and adds:
//   Global (3 eqs):
//     fggen ODE — total fission gas generated [mol]
//     fgrel ODE — total fission gas released [mol]
//     gpres AE  — ideal gas law: gpres = (mu0 + fgrel)*R/vt*1e-6 [MPa]
//                 where vt = (vgp + delta_annular_gap) / T_plenum
//                 (liquid-bond volume: sodium fills the gap; only the plenum
//                  plus the change from the initial annular gap volume are free)
//
//   Per-layer irradiation block (neq_irr = 1 + nf + nc):
//     bup_y = bup*8.64e10  ODE — burnup energy proxy
//     efs[0..nf-1]         ODEs — volumetric fuel swelling strains
//     ec[0..nc-1]          ODEs — cladding hoop creep strains (HT9::creepRate)
//
// y-vector layout:
//   [0]  fggen
//   [1]  fgrel
//   [2]  gpres
//   [3 + j*neq_j + 0..neq_th-1]                    thermal block
//   [3 + j*neq_j + neq_th]                          bup_y = bup*8.64e10
//   [3 + j*neq_j + neq_th+1..neq_th+nf]             efs[0..nf-1]
//   [3 + j*neq_j + neq_th+nf+1..neq_th+nf+nc]       ec[0..nc-1]
//   [3 + j*neq_j + neq_th+neq_irr..]                mechanics
// -----------------------------------------------------------------------
class FredMNaResiduals : public RodResiduals {
public:
    FredMNaResiduals(const FuelRodGeometry& geom,
                     UPuZr&                 fuel,
                     const HT9&             clad,
                     GapMaterial&           gap_mat,
                     double                 coolant_pressure_MPa);

    int neqGlobal() const { return 3; }
    int neqTotal()  const { return 3 + m_neq_j * m_geom.nz; }

    // MNa-specific boundary condition setters.
    void setInitialGasPressure(double p_MPa) {
        m_gpres0 = p_MPa;
        m_mu0 = p_MPa * 1.0e6 * m_geom.vgp / (R_GAS * T_REF);
    }
    void setBondPressure(double p_MPa) { setInitialGasPressure(p_MPa); } // legacy alias
    void setGasInventory(double mu0_mol) { m_mu0 = mu0_mol; }
    void setPlenumTemperature(TimeSeries fn) { m_plenumTFn = std::move(fn); }
    // Override: also update plenum temperature.
    void setInitialTemperature(double T0) override { m_T_init = T0; m_Tplenum = T0; }

    void   setElapsedTime(double t) { m_elapsed_time = t; }
    double elapsedTime()      const { return m_elapsed_time; }
    // Updated by FredMNaSolver::afterAcceptedStep from GRSIS cgrel output.
    void   setFgrelFromGrsis(double v) { m_fgrel_grsis = v; }
    double vfree0()           const { return m_vfree0; }   // plenum volume [m3]

    // FredMNaSolver::afterAcceptedStep computes irradiation physics (Zr
    // redistribution, GRSIS bubble swelling, irradiation conductivity
    // correction) on its OWN FredMNaRodState. Residual evaluation
    // (prepareLayer / computeThermalResiduals / computeIrradiationResiduals)
    // runs on this class's private m_state instead, which unpackState()
    // refreshes only from the IDA y-vector each Newton iteration — none of
    // the above physics is part of y, so it never crosses over on its own.
    // Call this once per accepted step (per layer) to push the auxiliary,
    // non-y-vector fields across; same one-step-lag pattern as
    // setFgrelFromGrsis, just generalized to the whole per-layer aux state.
    // Also copies flag + directional swelling (efsz/efsh/efsr) so
    // FredMNaStressStrain (which reads this class's own m_state, not the
    // solver's) sees the current gap-contact state (gap-behaviour-model
    // refactor, Decision D8).
    void syncAuxLayerState(int j, const FredMNaLayerState& src);

    // Gap state override (also updates layer flag string).
    void setGapOpen(int layer, bool open) override;

    // FRED-M-Na's own IDA residual entry point (Decision D2/D8 of the gap-
    // behaviour-model refactor). RodResiduals::computeResiduals (the shared
    // NVI loop) is non-virtual and calls the inherited StressStrain m_mech
    // directly, so it cannot be overridden without touching the platform
    // header (RodResiduals.hpp is untouched, Decision D1). Instead this
    // reimplements the same per-layer loop and calls FredMNaStressStrain
    // (three-state open/soft/clos BCs + directional swelling) in place of
    // m_mech. FredMNaSolver::getSolverResFn() registers idaResidualMNa
    // instead of RodResiduals::idaResidual.
    // \coderef{fred_m_na_residuals_mna_nvi}
    int computeResidualsMNa(double t, const double* y, const double* yp, double* r) const;
    static int idaResidualMNa(double t, N_Vector y, N_Vector yp, N_Vector r, void* user_data);

    // ---- One-step (backward-Euler, fixed dt) integrator support ----
    // FredMNaSolver::runOneStepLoop (Task: drop SUNDIALS IDA for FRED-M-Na,
    // replace with the legacy "backward-Euler, always accept" scheme) solves
    // each axial layer's local nonlinear block independently (layers are
    // coupled only through the 3 global rod-scalars fggen/fgrel/gpres, never
    // directly to each other), using a small per-layer dense Newton solve
    // instead of one global neqTotal()-sized Newton solve. These two methods
    // expose that block structure:
    //
    //   computeLayerResiduals — thermal + irradiation + mechanics residual
    //     for ONE layer only (O(1), not O(nz) like computeResidualsMNa),
    //     given that m_state's globals (gpres) and all OTHER layers are
    //     already current. Updates m_state.layers[j] from yj/ypj first.
    //   computeGlobalUpdate — the 3 global rows (fggen ODE, fgrel AE, gpres
    //     AE) are all explicit/direct-substitution given the CURRENT layer
    //     states (none of them implicitly depends on its own unknown), so
    //     they need no Newton solve at all; this returns the new values
    //     directly for the caller's fixed-point (Picard) outer sweep over
    //     globals + layers.
    void computeLayerResiduals(int j, double t, const double* yj, const double* ypj,
                                double* rj) const;
    void computeGlobalUpdate(double t, double dt, double fggen_old,
                              double& fggen_new, double& fgrel_new,
                              double& gpres_new) const;
    // Pushes the 3 global values (as solved by the caller's outer Picard
    // sweep, see computeGlobalUpdate) into m_state, so computeLayerResiduals'
    // reads of m_state.gpres (mechanics BC input) see the current value.
    // Without this, computeLayerResiduals never sees anything but the
    // globals' initial (zero) state, since it deliberately does NOT run a
    // full unpack() (that's the whole point of avoiding the O(nz) cost).
    void setGlobalState(double fggen, double fgrel, double gpres) const {
        m_state.fggen = fggen;
        m_state.fgrel = fgrel;
        m_state.gpres = gpres;
    }
    // Seeds m_state fully (all globals + all layers) from y/yp via unpack().
    // Call once before the one-step loop starts (cold start or restart) so
    // computeGlobalUpdate's very first call (before any computeLayerResiduals
    // call has had a chance to populate a layer) sees real geometry instead
    // of the zero-initialized default.
    void primeState(const double* y, const double* yp) const { unpack(y, yp, m_state); }
    // gpres from the CURRENT m_state (all layers, whatever their current
    // freshness). Used by solveStepBackwardEuler's joint layer+gpres Newton
    // (see newtonSolveLayer): gpres's own defining equation
    // (gpres = computeGasPressure(...)) has real cross-coupling to this
    // layer's own mechanics unknowns (efh/efz/eh/ez, via the fuel/clad
    // outer-radius and axial-strain terms in computeGasPressure), so it is
    // solved jointly with each layer rather than frozen as a constant
    // input, for a tighter within-step coupling than a separate direct-
    // substitution pass would give.
    double computeGasPressureCurrent() const { return computeGasPressure(m_state); }

    // Per-unknown "typical physical magnitude" for one layer's y-vector
    // slice (size neqPerLayer()) — used ONLY to floor the one-step
    // integrator's finite-difference Jacobian step (FredMNaDenseSolve.hpp),
    // not as a solution or convergence scale. This system mixes wildly
    // different physical units in a single dense residual block (e.g. qqv
    // ~1e9 W/m3 next to a strain ~1e-3), and several of the large-scale
    // unknowns (qqv, bup_y) legitimately start at exactly 0 (before power
    // ramp-up / at the first ever step) — a plain relative FD step
    // (eps*|current value|) silently underflows to a bitwise-zero
    // derivative in that case (the perturbation is smaller than one ULP of
    // the ~1e9-scale residual it feeds), which was found to make the very
    // first Newton solve of a cold-started run silently a no-op. Mirrors
    // fillIdArray's structure/offsets.
    void fillFdScale(double* scale) const;

    // Refreshes the subchannel coolant field (T_co / HTC per layer) for
    // time t. computeResidualsMNa calls this once per full residual
    // evaluation (before its per-layer loop); computeLayerResiduals is
    // scoped to a single layer and does NOT call it (that's the point of
    // avoiding the O(nz) cost), so the one-step integrator must call this
    // itself once per outer Picard sweep, or layerHTC()/the Robin outer-
    // clad-surface BC silently reads a stale/uninitialised coolant field.
    void updateSubchannelField(double t) const;

    // Pack / unpack between FredMNaRodState and the flat IDA vector.
    void pack  (const FredMNaRodState& state, double* y) const;
    void unpack(const double* y, const double* yp, FredMNaRodState& state) const;

    // id[] array for SUNDIALS IDASetId (1=ODE, 0=AE).
    void fillIdArray(double* id) const;

    // Gas pressure from ideal gas law with liquid-bond V/T: (mu0+fgrel)*R/vT*1e-6 [MPa].
    double computeGasPressure(const FredMNaRodState& state) const;

    // Initialise algebraic state variables to consistent starting values.
    void initAlgebraicState(FredMNaRodState& state) const;

protected:
    // ---- RodResiduals virtual hooks ----
    void   unpackState    (const double* y, const double* yp)  const override;
    int    globalOffset   ()                                    const override { return 3; }
    double globalPressure ()                                    const override { return m_state.gpres; }
    const AxialLayerState& layerState(int j)                   const override { return m_state.layers[j]; }
    void   prepareLayer   (int j)                              const override;

    void computeGlobalResiduals     (double t, const double* yp, double* r)                             const override;
    void computeThermalResiduals    (int j, double t, const double* yj, const double* ypj, double* rj)  const override;
    void computeIrradiationResiduals(int j, const double* ypj, double* rj)                              const override;
    void onGapClosed(int layer)                                                                               override;

private:
    UPuZr&             m_fuel;
    const HT9&         m_ht9;     // concrete cladding ref for HT9-specific methods
    const GapMaterial& m_gap_mat;

    double m_gpres0  = 0.1;
    double m_mu0     = 0.0;
    double m_vfree0  = 0.0;   // plenum volume = geom.vgp (liquid bond: gas confined to plenum)
    double m_h_total = 0.0;   // total initial active length = sum_j dz0[j]
    mutable double m_Tplenum = T_REF;
    double m_elapsed_time = 0.0;

    TimeSeries m_plenumTFn;

    double m_fgrel_grsis = 0.0;  // fgrel set by GRSIS; AE residual drives y[1] here

    mutable FredMNaRodState m_state;

    // Per-layer slice of unpack() — sets one FredMNaLayerState from its
    // yj/ypj slice of the flat y/yp vector. Factored out so it can also be
    // called for a single layer (computeLayerResiduals) without re-unpacking
    // every other layer (see one-step integrator support above).
    void unpackLayer(int j, const double* yj, const double* ypj, FredMNaLayerState& s) const;

    // FRED-M-Na's own three-state mechanics (Decision D2/D8); used only by
    // computeResidualsMNa, in place of the inherited (two-state) m_mech.
    FredMNaStressStrain         m_mech_mna;
    mutable std::vector<double> m_r_mech_mna;

    // Intra-layer y-vector offsets (relative to layer start).
    int offBup()      const { return m_neq_th; }
    int offEfs(int i) const { return m_neq_th + 1 + i; }
    int offEc (int i) const { return m_neq_th + 1 + m_geom.nf + i; }
    int offMech()     const { return m_neq_th + m_neq_irr; }

    // Mechanics variable offsets within the mechanics block.
    int offEft  (int i) const { return offMech() + i; }
    int offEfh  (int i) const { return offMech() + m_geom.nf + i; }
    int offEfr  (int i) const { return offMech() + 2*m_geom.nf + i; }
    int offEfz  ()      const { return offMech() + 3*m_geom.nf; }
    int offSigfh(int i) const { return offMech() + 3*m_geom.nf+1 + i; }
    int offSigfr(int i) const { return offMech() + 4*m_geom.nf+1 + i; }
    int offSigfz(int i) const { return offMech() + 5*m_geom.nf+1 + i; }
    int offEt   (int i) const { return offMech() + 6*m_geom.nf+1 + i; }
    int offEh   (int i) const { return offMech() + 6*m_geom.nf+m_geom.nc+1 + i; }
    int offEr   (int i) const { return offMech() + 6*m_geom.nf+2*m_geom.nc+1 + i; }
    int offEz   ()      const { return offMech() + 6*m_geom.nf+3*m_geom.nc+1; }
    int offSigh (int i) const { return offMech() + 6*m_geom.nf+3*m_geom.nc+2 + i; }
    int offSigr (int i) const { return offMech() + 6*m_geom.nf+4*m_geom.nc+2 + i; }
    int offSigz (int i) const { return offMech() + 6*m_geom.nf+5*m_geom.nc+2 + i; }
};

} // namespace fred
