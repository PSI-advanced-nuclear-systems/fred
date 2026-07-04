#pragma once
#include "../../platform/RodResiduals.hpp"
#include "../../platform/GapPressureModel.hpp"
#include "platform/GapMaterial.hpp"

namespace fred {

class FuelPelletMaterial;
class CladdingMaterial;

// \coderef{fred_rod_residuals}
// -----------------------------------------------------------------------
// FredRodResiduals — DAE residual assembler for FRED-ROD.
//
// Inherits shared infrastructure from RodResiduals and adds:
//   • 1 global equation: gas pressure AE (delegated to GapPressureModel)
//   • Per-layer: thermal block + mechanics block (no irradiation)
//
// \coderef{fred_rod_yvec_layout}
// y-vector layout:
//   y[0]                                    = gpres  (global AE, id=0)
//   y[1 + j*neq_j .. 1+(j+1)*neq_j - 1]    = layer j (thermal + mechanics)
//
// Within each layer (neq_j = neq_th + neq_mech = 6 + 7*nf + 7*nc):
//   Thermal block  (neq_th = 4 + nf + nc):  qqv, gap, pfc, hgap, T_fuel, T_clad
//   Mechanics block (neq_mech = 6nf+6nc+2): eft,efh,efr,efz,sigfh,sigfr,sigfz,
//                                             et,eh,er,ez,sigh,sigr,sigz
// -----------------------------------------------------------------------
class FredRodResiduals : public RodResiduals {
public:
    FredRodResiduals(const FuelRodGeometry&    geom,
                     const FuelPelletMaterial& fuel,
                     const CladdingMaterial&   clad,
                     const GapMaterial&        gap_mat,
                     double coolant_pressure_MPa);

    // Total equations: 1 global (gpres) + nz per-layer blocks.
    int neqTotal() const { return 1 + m_neq_j * m_geom.nz; }

    // Inject the gap pressure model (non-owning; caller retains ownership).
    void setGapPressureModel(GapPressureModel& model) { m_pressure_model = &model; }

    // Set algebraic state variables to physically consistent initial values.
    void initAlgebraicState(RodState& state,
                            const FuelPelletMaterial& fuel,
                            const CladdingMaterial& clad) const;

    // Newton solve for FVM-consistent mechanical initial conditions.
    void solveMechanicalIC(RodState& state, double t, double rtol, double atol) const;

    // Pack / unpack between RodState and the SUNDIALS flat vector.
    void pack  (const RodState& state, double* y) const;
    void unpack(const double* y, const double* yp, RodState& state) const;

    // Called at gap REOPENING: updates gap manager and layer state.
    void setGapOpen(int layer, bool open) override;

protected:
    // ---- RodResiduals virtual hooks ----
    void unpackState(const double* y, const double* yp) const override;
    int  globalOffset()                                 const override { return 1; }
    double globalPressure()                             const override { return m_state.gpres; }
    const AxialLayerState& layerState(int j)            const override { return m_state.layers[j]; }

    // \coderef{fred_rod_thermal_hook}
    void computeGlobalResiduals  (double t, const double* yp, double* r)                             const override;
    void computeThermalResiduals (int j, double t, const double* yj, const double* ypj, double* rj)  const override;
    void freezeThermalResiduals  (int j, const double* yj, double* rj)                               const override;
    void onGapClosed             (int layer)                                                                 override;

private:
    const GapMaterial&  m_gap_mat;
    GapPressureModel*   m_pressure_model = nullptr; // non-owning

    mutable RodState    m_state;

    // y-variable offsets within a layer block (0-based, from layer start).
    int offEft()   const { return m_neq_th; }
    int offEfh()   const { return offEft()   + m_geom.nf; }
    int offEfr()   const { return offEfh()   + m_geom.nf; }
    int offEfz()   const { return offEfr()   + m_geom.nf; }
    int offSigfh() const { return offEfz()   + 1; }
    int offSigfr() const { return offSigfh() + m_geom.nf; }
    int offSigfz() const { return offSigfr() + m_geom.nf; }
    int offEt()    const { return offSigfz() + m_geom.nf; }
    int offEh()    const { return offEt()    + m_geom.nc; }
    int offEr()    const { return offEh()    + m_geom.nc; }
    int offEz()    const { return offEr()    + m_geom.nc; }
    int offSigh()  const { return offEz()    + 1; }
    int offSigr()  const { return offSigh()  + m_geom.nc; }
    int offSigz()  const { return offSigr()  + m_geom.nc; }
};

} // namespace fred
