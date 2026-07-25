#pragma once
#include "platform/CladdingMaterial.hpp"

namespace fred {

// FRED-M-Na cladding material interface: extends the shared CladdingMaterial
// base with the irradiation-response physics metallic-fuel SFR cladding
// needs (creep, void swelling, lanthanide-driven wastage) that FRED-ROD /
// FRED-OX's generic interface has no room for — those effects need extra
// rate-law inputs (qqv, elapsed time, dose, flux) that plain (T, sigma)
// can't carry.
//
// Model SELECTION is solver-level (CladSwellingModel / CladWastageModel,
// set via FredMNaSolver::set_clad_swelling_model / set_cladding_wastage_model):
// FRED-M-Na fixes which numerical scheme evaluates each named effect
// (Hofmann's closed-form fluence function, SAS4A's dose-gated piecewise
// rate, the precipitation-kinetics front law, the lanthanide-diffusion
// PDE — see FredMNaCladSwelling.hpp / FredMNaCladWastage.hpp for the
// generic scheme code). A material subclass supplies only the
// correlation values a specific alloy plugs into those fixed schemes —
// it does not get to redefine the schemes themselves. That is a
// deliberate scope boundary: a cladding alloy with a genuinely different
// physical mechanism (not just different numbers) belongs in a new
// application built on top of FRED-M-Na, not as ever-more-generic
// scaffolding bolted onto this interface.
//
// All new virtuals default to a physically-inert value (0.0 — "this
// effect is not modelled for this material"), so a minimal subclass
// compiles and runs with irradiation effects opted out, matching
// FredMNaSolver's own default-off philosophy for these features
// (set_enable_clad_swelling / set_enable_clad_wastage default false).
class FredMNaCladdingMaterial : public CladdingMaterial {
public:
    // --- Failure margins (FredMNaFailure.hpp: burst / melt criteria) ---
    virtual double yieldStress(double T_K) const = 0;
    virtual double burstStress(double T_K, double ssy0_Pa) const = 0;

    // --- Thermal + irradiation creep rate [1/s] (FredMNaResiduals.cpp,
    //     clad hoop creep ODE) ---
    //   sig_Pa : effective hoop stress [Pa]
    //   qqv    : volumetric power density [W/m3] (neutron-flux proxy)
    //   time_s : elapsed time [s] (tertiary-creep term)
    virtual double creepRate(double T_K, double sig_Pa,
                              double qqv, double time_s) const {
        (void)T_K; (void)sig_Pa; (void)qqv; (void)time_s;
        return 0.0;
    }

    // --- Void swelling: linear (volumetric/3) strain contribution ---
    // Hofmann model: closed-form function of cumulative fast fluence,
    // evaluated directly (no external integration needed).
    //   neuflue : fast fluence [1e22 n/cm2]
    // Returns the volumetric strain [-] (caller applies the /3 split).
    virtual double voidSwelling(double neuflue, double T_K) const {
        (void)neuflue; (void)T_K;
        return 0.0;
    }
    // SAS4A model: piecewise-in-T strain RATE [1/s]. The solver applies
    // the dose-onset gate (CladSwellingParams::dose_threshold_dpa)
    // generically before calling this — implementations only need to
    // return their own T-dependent rate curve for the gated regime.
    //   flux_1e22 : instantaneous fast flux [1e22 n/cm2/s]
    //   dconv     : fluence->dose conversion factor [dpa per 1e22 n/cm2]
    //   dose_dpa  : cumulative cladding dose [dpa] (already >= threshold)
    virtual double voidSwellingSAS4ARate(double T_K, double flux_1e22,
                                          double dconv, double dose_dpa) const {
        (void)T_K; (void)flux_1e22; (void)dconv; (void)dose_dpa;
        return 0.0;
    }

    // --- Cladding wastage correlations (FredMNaCladWastage.hpp) ---
    // Precipitation-kinetics (PK) model:
    virtual double wastageDiffusionPK(double T_K) const { (void)T_K; return 0.0; }  // D(T) [m2/s]
    virtual double wastageSolubilityFuel() const { return 0.0; }  // c_sol,fuel [-]
    virtual double wastageSolubilityClad() const { return 0.0; }  // c_sol,clad [-]
    // Lanthanide-tracking (LaTracking) model:
    virtual double wastageDiffusionLaTracking(double T_K) const { (void)T_K; return 0.0; } // D_La(T) [m2/s]
    virtual double wastageUptakeCapacity() const { return 0.0; }  // u_sol [#/m3]
};

} // namespace fred
