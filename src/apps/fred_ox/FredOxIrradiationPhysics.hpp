#pragma once
// FredOxIrradiationPhysics.hpp
//
// Self-contained irradiation physics functions for the FRED-OX application.
// Kept in a separate file so the irradiation model can be inspected and
// modified independently of the main solver.
//
// All functions are free functions (no class state) matching the formulas
// in legacy FRED.f90 exactly.  See session_state.md for equation references.
//
// Units convention (same as legacy FRED globals.f90):
//   bup      [MWd/kgU]     — local fuel burnup
//   qqv      [W/m3]        — local volumetric power density
//   rof0     [kg/m3]       — initial fuel density
//   T        [K]           — temperature
//   rate     [fiss/m3/s]   — fission rate density

#include <cmath>
#include <algorithm>
#include <vector>

namespace fred {

// Fission energy per fission [J] = 210 MeV * 1.60214e-13 J/MeV
constexpr double ENERGY_PER_FISSION = 210.0 * 1.60214e-13; // J/fiss

// Avogadro number
constexpr double AVOGADRO_OX = 6.0247e23; // 1/mol

// Burnup scaling: 1 MWd/kgU = 8.64e10 J/kgU
// The y-vector stores bup * BUP_SCALE to keep ODE in W/kgU units
// (so the residual dbup_y/dt = qqv/rof0 is dimensionally consistent).
constexpr double BUP_SCALE = 8.64e10; // J per MWd/kgU

// -------------------------------------------------------------------------
// Fission rate density  [fissions / (m3·s)]
// -------------------------------------------------------------------------
inline double fissionRate(double qqv_W_m3) {
    return qqv_W_m3 / ENERGY_PER_FISSION;
}

// -------------------------------------------------------------------------
// Waltar-Reynolds fission gas release fraction for a single axial layer.
//
// Returns the dimensionless fraction [0..1] of generated fission gas
// that has been released to the free volume for a given layer.
//
//   bup      : local burnup [MWd/kgU]
//   T_nodes  : temperature at each fuel radial node [K], size nf
//   area0    : cross-sectional area of each fuel ring [m2], size nf
//              (area0[i] = π * drf0 * rad0[i] for interior rings)
//   nf       : number of fuel radial nodes
//   pucont   : Pu mole fraction [-]  (used for melting temperature)
//   ql_kWm   : linear power density [W/m] for this layer (used for fgrU)
//
// Reference: Waltar & Reynolds (1981), fast breeder reactor, p.198.
// -------------------------------------------------------------------------
inline double fgReleaseRate(double bup, const double* T_nodes,
                             const double* area0, int nf,
                             double pucont, double ql_Wm)
{
    if (bup <= 0.0) return 0.0;

    // Melting temperature for MOX [K]
    const double tM = 3120.0 - 388.1 * pucont - 30.4 * pucont * pucont;

    // --- Release fraction in restructured zone ---
    // Valid for any burnup > 0.
    const double a_R = 4.7 / bup * (1.0 - std::exp(-bup / 5.9));
    const double fgrR = std::max(0.0, 1.0 - a_R);

    // --- Release fraction in unrestructured zone ---
    double fgrU = 0.0;
    if (bup > 3.5) {
        const double ql_kWm = ql_Wm / 1000.0;
        double a_U = 25.6 / (bup - 3.5) * (1.0 - std::exp(-(bup / 3.5 - 1.0)))
                     * std::exp(-0.0125 * ql_kWm);
        if (bup >= 49.2)
            a_U *= std::exp(-0.3 * (bup - 49.2));
        fgrU = std::max(0.0, 1.0 - a_U);
    }

    // --- Fractional areas for each zone ---
    // T > tM-50 K  → molten zone (fgr = 1)
    // T > 1773 K (1500°C) and T < tM-50 → restructured zone
    // else → unrestructured
    double aM = 0.0, aR = 0.0, aU = 0.0, aTotal = 0.0;
    for (int i = 0; i < nf; ++i) {
        aTotal += area0[i];
        if (T_nodes[i] >= tM - 50.0)
            aM += area0[i];
        else if (T_nodes[i] - 273.15 >= 1500.0)
            aR += area0[i];
        else
            aU += area0[i];
    }

    if (aTotal <= 0.0) return 0.0;
    return (aM * 1.0 + aR * fgrR + aU * fgrU) / aTotal;
}

// -------------------------------------------------------------------------
// MATPRO fuel volumetric swelling rate  [1/s]
//
// Returns the volumetric swelling strain rate (d(ε_s)/dt) in [1/s].
// The caller multiplies by a user-supplied multiplier (default 1.0)
//
//   T        : local fuel temperature [K]
//   bup      : local burnup [MWd/kgU]
//   qqv      : local power density [W/m3]
//   rof0     : initial fuel density [kg/m3]
//
// Reference: MATPRO-09 (NUREG/CR-0479), as implemented in FRED.f90 fswel.
// Gas FP swelling + solid FP swelling.
// -------------------------------------------------------------------------
inline double fuelSwellingRate(double T_K, double bup_MWdkgU,
                                double qqv_W_m3, double rof0_kg_m3)
{
    const double rate = qqv_W_m3 / ENERGY_PER_FISSION; // fiss/(m3·s)

    double fswel = 0.0;

    // Gaseous fission product swelling (MATPRO model)
    if (T_K < 2800.0) {
        const double dT = 2800.0 - T_K;
        // Total energy deposited in fuel [J/m3]
        const double E_Jm3 = bup_MWdkgU * BUP_SCALE * rof0_kg_m3;
        // Total fissions per unit volume [fiss/m3]
        const double B = E_Jm3 / ENERGY_PER_FISSION;
        fswel += 8.8e-56 * std::pow(dT, 11.73)
                 * std::exp(-0.0162 * dT)
                 * std::exp(-8.0e-27 * B)
                 * rate;
    }

    // Solid fission product swelling (MATPRO)
    fswel += 2.5e-29 * rate;

    return fswel;
}

// -------------------------------------------------------------------------
// Fuel relocation fraction (FRAPCON model)
//
//   ql_Wm   : linear power [W/m]
//   bup     : burnup [MWd/kgU]
//
// Returns 0 when FUEL_RELOC is not enabled.
// -------------------------------------------------------------------------
inline double fuelRelocFraction(double ql_Wm, double bup) {
    const double ql_kWm = ql_Wm / 1000.0;
    return 0.25 + std::max(5.0, std::min(10.0, ql_kWm / 4.0))
           * (std::min(1.0, bup / 5.0) + 1.0) / 100.0;
}

} // namespace fred
