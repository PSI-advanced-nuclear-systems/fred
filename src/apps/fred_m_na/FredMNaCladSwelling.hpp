#pragma once
// Cladding void-swelling models for FRED-M-Na.
//
// Two selectable models (FredMNaSolver::setCladSwellingModel); the
// numerical scheme for each is fixed here (application-level), while the
// per-material correlation values come from the FredMNaCladdingMaterial
// passed to the solver (see FredMNaCladdingMaterial.hpp: voidSwelling(),
// voidSwellingSAS4ARate()):
//
//   Hofmann (default) — Cswel.for icswel=1, Hofmann (1985)-shaped
//     void-swelling correlation.  Strain is evaluated directly from the
//     cumulative fast fluence (no time integration needed), per radial
//     clad node at that node's own temperature.
//
//   SAS4A — Cswel.for icswel=2 (added Timpano/Daniele, 2024-08-06):
//     piecewise-in-temperature swelling RATE, gated on cladding dose
//     >= dose_threshold_dpa (default 100 dpa, applied here — see the
//     dose-gate call site in FredMNaSolver.cpp), evaluated once at the
//     cladding mid-wall temperature and integrated (explicit Euler) into a
//     strain, applied uniformly across the clad wall (legacy: "use
//     midwall temperature if model 2 is selected").
//
// Both models consume a common SAS-derived fast-neutron fluence/dose
// bookkeeping (Baseir.for lines ~155-162):
//   flux [1e22 n/cm2/s] = ql * fltpow * 1e-26,   ql = qqv * pi*(rfo0^2-rfi0^2)
//   neuflue2 (fluence, 1e22 n/cm2)  += flux * dt
//   dose (dpa)                      += flux * dconv * dt
// NOTE: fltpow (fast-flux-to-linear-power ratio) is a case-specific
// neutronics calibration constant (Serpent-derived, legacy comment: "fast
// flux over linear power ratio from Serpent for Inner Fuel Average
// conditions"), not a universal physical constant.  The default here is
// Timpano's Inner-Fuel-Average value; re-derive it (or at least sanity
// check it) for other pin/core positions via
// FredMNaSolver::setFastFluxCalibration().

namespace fred {

enum class CladSwellingModel {
    Hofmann = 1,   // Cswel.for icswel=1
    SAS4A   = 2    // Cswel.for icswel=2
};

struct CladSwellingParams {
    double fltpow = 6.13209e14;   // fast flux / linear power [Timpano IF-AVG Serpent fit]
    double dconv  = 4.0;          // fluence[1e22 n/cm2] -> dose[dpa], SAS User's Manual (Superphenix)
    double dose_threshold_dpa = 100.0;  // SAS4A swelling onset [dpa] (scheme-level gate,
                                        // applied before calling into the material)
};

} // namespace fred
