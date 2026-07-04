#pragma once
#include "FuelRodGeometry.hpp"
#include "FuelRodState.hpp"
#include "GapMaterial.hpp"
#include "Constants.hpp"
#include <functional>
#include <algorithm>

namespace fred {

// -----------------------------------------------------------------------
// computeGasBondVT — temperature-weighted free volume for a gas-bond pin.
//
// Sums contributions from:
//   plenum:       V_gp / T_plenum
//   annular gap:  sum_j pi * dz_j * (rci_j^2 - rfo_j^2) / T_gap_j
//   central void: sum_j pi * dz_j * rfi_j^2 / T_ctr_j   (if rfi0 > 0)
//
// LayerContainer elements must expose T[], rci, rfo, rfi (AxialLayerState or
// any struct that derives from it, e.g. FredOxLayerState).
// -----------------------------------------------------------------------
template<typename LayerContainer>
double computeGasBondVT(const FuelRodGeometry& geom,
                         const LayerContainer& layers,
                         double T_plenum)
{
    const int nf = geom.nf;
    double vt = geom.vgp / std::max(T_plenum, 1.0);
    for (int j = 0; j < geom.nz; ++j) {
        const auto& s = layers[j];
        double T_gap = 0.5 * (s.T[nf-1] + s.T[nf]);
        if (T_gap < 1.0) T_gap = 1.0;
        vt += PI * geom.dz0[j] * (s.rci*s.rci - s.rfo*s.rfo) / T_gap;
        if (geom.rfi0 > 0.0) {
            double T_ctr = s.T[0];
            if (T_ctr < 1.0) T_ctr = 1.0;
            vt += PI * geom.dz0[j] * s.rfi*s.rfi / T_ctr;
        }
    }
    return vt;
}

// -----------------------------------------------------------------------
// computeLiquidBondVT — temperature-weighted free volume for a liquid-bond pin.
//
// For sodium-bonded (liquid-bond) pins the fill gas is confined to the plenum.
// The effective gas volume changes as the annular gap opens or closes:
//   V_T = (V_plenum + delta_annular_vol) / T_plenum
// where
//   delta_annular_vol = sum_j pi*dz_j*(rci_j^2 - rfo_j^2)  — current gap vol
//                     - pi*(rci0^2 - rfo0^2)*h_total         — initial gap vol
//
// h_total must equal sum_j dz0[j] (total initial active length).
//
// LayerContainer elements must expose rci, rfo (AxialLayerState or any
// struct that derives from it, e.g. FredMNaLayerState).
// -----------------------------------------------------------------------
template<typename LayerContainer>
double computeLiquidBondVT(const FuelRodGeometry& geom,
                            const LayerContainer& layers,
                            double T_plenum, double h_total)
{
    const double gap_area0 = PI * (geom.rci0*geom.rci0 - geom.rfo0*geom.rfo0);
    double vt = geom.vgp;
    for (int j = 0; j < geom.nz; ++j) {
        const auto& s = layers[j];
        vt += PI * geom.dz0[j] * (s.rci*s.rci - s.rfo*s.rfo);
    }
    vt -= gap_area0 * h_total;
    return vt / std::max(T_plenum, 1.0);
}

// -----------------------------------------------------------------------
// GapPressureModel — residual for the global gas pressure algebraic equation.
//
// Gas bond (selected by the gas-bond constructor):
//   Ideal gas law for the fill-gas inventory.
//   Handles gas-bond and liquid-bond pins via GapMaterial::isGasBond():
//     Gas bond:    v_T = computeGasBondVT(geom, layers, T_plenum)
//     Liquid bond: v_T = (V_plenum + delta_annular_vol) / T_plenum
//   Residual: r = gpres - mu0 * R / v_T * 1e-6
// -----------------------------------------------------------------------
// \coderef{gap_pressure_model}
class GapPressureModel {
public:
    // Gas bond constructor: calculates mu0 from the as-fabricated free volume.
    GapPressureModel(double gpres0_MPa, const FuelRodGeometry& geom, const GapMaterial& gapmat);
    virtual ~GapPressureModel() = default;

    // Change the initial fill-gas pressure and recompute the gas inventory.
    // Only valid when constructed with the gas-bond constructor.
    void setInitialGasPressure(double gpres0_MPa);

    virtual double residual(double t, double gpres,
                            const RodState& state,
                            const FuelRodGeometry& geom) const;
    virtual double initialPressure() const;
    // Notify the model of the initial temperature [K] (sets plenum temperature).
    virtual void setInitialTemperature(double T0);
    // Optional time-varying plenum temperature; overrides setInitialTemperature
    // after t=0.  Call before run().
    void setPlenumTemperature(std::function<double(double)> fn);

private:
    // Gas bond fields
    double m_gpres0          = 0.1;
    double m_mu0             = 0.0;
    mutable double m_Tplenum = T_REF;   // mutable: updated lazily in residual()
    double m_h_total         = 0.0;
    const FuelRodGeometry* m_geom   = nullptr;
    const GapMaterial*     m_gapmat = nullptr;
    std::function<double(double)> m_plenumTFn;
};

} // namespace fred
