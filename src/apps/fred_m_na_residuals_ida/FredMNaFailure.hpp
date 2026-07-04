#pragma once
// FredMNaFailure.hpp — Failure criteria for FRED-M-Na (U-Pu-Zr / HT-9 pin)
//
// Ported from Cfail.for, Csigb.for, Mmelt.for (Timpano/EPFL FRED_M_OCT24.SRC).
//
// Criterion 1 — Cladding burst:
//   fcrit_burst = sigma_eff / burst_stress
//   burst_stress from HT9::burstStress(T, ssy0)
//
// Criterion 4 — Fuel melt margin:
//   fcrit_melt = T_fuel / T_solidus
//   T_solidus from UPuZr::solidusTemperature(zr_wf)
//
// Both criteria return values in [0, 1]; fcrit >= 1 signals failure.

#include "apps/fred_m_na/claddingmaterial/HT9.hpp"
#include "apps/fred_m_na/fuelpelletmaterial/UPuZr.hpp"

namespace fred {

// -------------------------------------------------------------------------
// FredMNaFailureState — tracked failure metrics for one pin at one time step
// -------------------------------------------------------------------------
struct FredMNaFailureState {
    // Per-layer values — size nz
    std::vector<double> fcrit_burst; // crit-1: sigma_eff / sigma_burst [-]
    std::vector<double> fcrit_melt;  // crit-4: T_fuel_peak / T_solidus [-]

    // Scalar summary flags
    bool   failed_burst = false;
    bool   failed_melt  = false;
    int    burst_layer  = -1;   // layer index where burst criterion first reached 1
    int    melt_layer   = -1;

    void resize(int nz) {
        fcrit_burst.assign(nz, 0.0);
        fcrit_melt .assign(nz, 0.0);
    }
};

// -------------------------------------------------------------------------
// computeFailureCriteria — evaluate burst and melt margins for all layers.
//
//   clad      : HT9 object (for yieldStress and burstStress)
//   fuel      : UPuZr object (for solidusTemperature)
//   nz        : number of axial layers
//   nf        : number of fuel radial nodes per layer
//   T_fuel[]  : per-layer peak fuel temperature [K], size nz (T at node 0)
//   T_clave[] : per-layer average cladding temperature [K], size nz
//   sigma_eff[]: per-layer effective cladding stress (volume-weighted) [Pa], size nz
//   zr_wf[]   : per-layer Zr weight fraction [-], size nz (from node-average)
//
// Updates `out` in-place.
// -------------------------------------------------------------------------
inline void computeFailureCriteria(const HT9&  clad,
                                    const UPuZr& fuel,
                                    int nz,
                                    const double* T_fuel,
                                    const double* T_clave,
                                    const double* sigma_eff_Pa,
                                    const double* zr_wf,
                                    FredMNaFailureState& out)
{
    out.resize(nz);
    out.failed_burst = false;
    out.failed_melt  = false;
    out.burst_layer  = -1;
    out.melt_layer   = -1;

    for (int j = 0; j < nz; ++j) {
        // Criterion 1: burst
        const double T_c   = T_clave[j];
        const double ssy0  = clad.yieldStress(T_c);
        const double sigb  = clad.burstStress(T_c, ssy0);
        out.fcrit_burst[j] = (sigb > 0.0) ? std::min(sigma_eff_Pa[j] / sigb, 1.0) : 0.0;

        if (out.fcrit_burst[j] >= 1.0 && !out.failed_burst) {
            out.failed_burst = true;
            out.burst_layer  = j;
        }

        // Criterion 4: fuel melt
        const double T_sol = fuel.solidusTemperature(zr_wf[j]);
        out.fcrit_melt[j]  = (T_sol > 0.0) ? std::min(T_fuel[j] / T_sol, 1.0) : 0.0;

        if (out.fcrit_melt[j] >= 1.0 && !out.failed_melt) {
            out.failed_melt = true;
            out.melt_layer  = j;
        }
    }
}

} // namespace fred
