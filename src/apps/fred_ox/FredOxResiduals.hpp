#pragma once
#include "gapmaterial/FredOxGapMaterial.hpp"
#include "FredOxState.hpp"
#include "platform/RodResiduals.hpp"
#include "platform/HeatConduction.hpp"
#include "fuelpelletmaterial/FredOxMOX.hpp"
#include <memory>
#include <vector>

namespace fred {

class CladdingMaterial;

// \coderef{fred_ox_residuals}
// -----------------------------------------------------------------------
// FredOxResiduals — DAE residual assembler for FRED-OX.
//
// Inherits shared infrastructure from RodResiduals and adds:
//   Global (3 eqs):
//     fggen ODE — total fission gas generated [mol/s]
//     fgrel ODE — total fission gas released [mol/s]
//     gpres AE  — gas pressure from ideal gas law with fission gas
//
//   Per layer (neq_irr = 1 + nf + nc extra equations):
//     bup     ODE  — burnup:          d(bup*8.64e10)/dt = qqv/rof0
//     efs[nf] ODEs — fuel swelling:   d(efs[i])/dt = fuelSwellingRate(...)
//     ec[nc]  ODEs — clad hoop creep: d(ec[i])/dt  = creepRate(T[nf+i], sigh[i])
//     Mechanics: same as FRED-ROD but Hooke's law includes swelling term.
//
// y-vector layout:
//   [0]  fggen
//   [1]  fgrel
//   [2]  gpres
//   [3 + j*neq_j + 0..neq_th-1]           thermal block
//   [3 + j*neq_j + neq_th]                        bup_y = bup*8.64e10
//   [3 + j*neq_j + neq_th+1..neq_th+nf]           efs[0..nf-1]
//   [3 + j*neq_j + neq_th+nf+1..neq_th+nf+nc]     ec[0..nc-1]
//   [3 + j*neq_j + neq_th+neq_irr..]               mechanics
// -----------------------------------------------------------------------
class FredOxResiduals : public RodResiduals {
public:
    FredOxResiduals(const FuelRodGeometry& geom,
                    FredOxMOX&             mox,
                    const CladdingMaterial& clad,
                    FredOxGapMaterial&      gap_ox,
                    double coolant_pressure_MPa,
                    double rof0_kg_m3,
                    double pu_content,
                    double sto0,
                    double swelling_multiplier = 1.0);

    int neqGlobal() const { return 3; }
    int neqTotal()  const { return 3 + m_neq_j * m_geom.nz; }

    // FRED-OX specific boundary condition setters.
    void setPlenumTemperature(TimeSeries fn)  { m_plenumTFn = std::move(fn); }
    void setInitialGasPressure(double p)      { m_gpres0 = p; m_gpres = p; }
    void setGasInventory(double mu0_mol)      { m_mu0 = mu0_mol; }
    void setSwellingMultiplier(double fswm)   { m_fswelmlt = fswm; }
    void setPuContent(double pu)              { m_pu = pu; }
    void setSto0(double sto)                  { m_sto0 = sto; }
    // Also update plenum temperature when the initial temperature is set.
    void setInitialTemperature(double T0) override { m_T_init = T0; m_Tplenum = T0; }

    // Gap state — setGapClosed needs type-specific unpack.
    void setGapOpen(int layer, bool open) override;

    // Pack / unpack between FredOxRodState and the SUNDIALS flat vector.
    void pack  (const FredOxRodState& state, double* y) const;
    void unpack(const double* y, const double* yp, FredOxRodState& state) const;

    // id[] array for SUNDIALS IDASetId (1=ODE, 0=AE).
    void fillIdArray(double* id) const;

    // Gas pressure: gpres = (mu0 + fgrel) * R / v_T * 1e-6 [MPa]
    double computeGasPressure(const FredOxRodState& state) const;

    // As-fabricated free volume [m3].
    double vfree0() const { return m_vfree0; }

    // Initialise algebraic state variables to physically consistent values.
    void initAlgebraicState(FredOxRodState& state) const;

    // Newton solve for mechanical initial conditions.
    void solveMechanicalIC(FredOxRodState& state, double t,
                           double rtol, double atol) const;

protected:
    // ---- RodResiduals virtual hooks ----
    void unpackState(const double* y, const double* yp)  const override;
    int  globalOffset()                                   const override { return 3; }
    double globalPressure()                               const override { return m_state.gpres; }
    const AxialLayerState& layerState(int j)              const override { return m_state.layers[j]; }
    void prepareLayer(int j)                              const override;

    // \coderef{fred_ox_global_residuals}
    void computeGlobalResiduals     (double t, const double* yp, double* r)                             const override;
    void computeThermalResiduals    (int j, double t, const double* yj, const double* ypj, double* rj)  const override;
    // \coderef{fred_ox_irradiation_hook}
    void computeIrradiationResiduals(int j, const double* ypj, double* rj)                              const override;
    void onGapClosed                (int layer)                                                                override;

    // Thread-safety: FredOxMOX/FredOxGapMaterial are mutated (setBurnup,
    // setLinearPower, ...) then read back within computeThermalResiduals,
    // so each OpenMP worker thread needs its own clone (see RodResiduals's
    // prepareThreadLocalState doc comment).
    void prepareThreadLocalState(int nthreads) const override;

private:
    // OX-specific material references (non-const — setters called during residual assembly).
    FredOxMOX&         m_mox;
    FredOxGapMaterial& m_gap_ox;

    // OX-specific fuel parameters.
    double m_rof0;      // initial fuel density [kg/m3]
    double m_pu;        // Pu content [-]
    double m_sto0;      // initial stoichiometry
    double m_fswelmlt;  // swelling multiplier

    // OX-specific gas state.
    double m_gpres0 = 0.1;
    double m_gpres  = 0.1;
    double m_mu0    = 0.0;
    double m_vfree0 = 0.0;
    mutable double m_Tplenum = T_REF;

    // OX-specific per-layer boundary conditions.
    TimeSeries m_plenumTFn;

    mutable FredOxRodState m_state;

    // Per-OpenMP-thread clones of the mutable material objects (see
    // prepareThreadLocalState). Index [0] is only ever used for threads==1.
    mutable std::vector<FredOxMOX>                  m_mox_tls;
    mutable std::vector<FredOxGapMaterial>           m_gap_ox_tls;
    mutable std::vector<std::unique_ptr<HeatConduction>> m_heat_tls;

    // Intra-layer offsets (relative to layer start in y).
    // Irradiation block: [bup | efs[0..nf-1] | ec[0..nc-1]]
    int offBup()       const { return m_neq_th; }
    int offEfs(int i)  const { return m_neq_th + 1 + i; }
    int offEc (int i)  const { return m_neq_th + 1 + m_geom.nf + i; }
    int offMech()      const { return m_neq_th + m_neq_irr; }

    // Mechanics offsets within the mech block (same relative layout as FRED-ROD).
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
