#pragma once
// FredMNaIrradiationPhysics.hpp
//
// Free physics functions for the FRED-M-Na application (U-Pu-Zr metallic fuel
// in sodium-cooled fast reactor).  All correlations are ported directly from
// the Fortran FRED_M_OCT24.SRC reference code without modification to the
// equations.  See doc/physics_fred_m_na.md for equation references.
//
// Sources:
//   Karahan, A. (2009). Modelling of thermo-mechanical and irradiation
//     behavior of metallic and oxide fuels for sodium fast reactors.
//   Hofman (1985). Metallic fuels handbook.
//   Ogata (1999). ALFUS irradiation behavior code.
//   SAS4A Users Manual.
//   Baseir.for / Sinf.for / Zrdist.for / FGR.for / Clanth.for (Timpano/EPFL)

#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

namespace fred {

// -------------------------------------------------------------------------
// Constants
// -------------------------------------------------------------------------
constexpr double MNA_AVOGADRO = 6.0247e23;
constexpr double MNA_R_GAS    = 8.314e-3; // kJ/(mol·K) — used in Zrdist diffusion
constexpr double MNA_R_CAL    = 1.98722;  // cal/(mol·K) — used in creep models
constexpr double MNA_PI       = 3.14159265358979323846;

// -------------------------------------------------------------------------
// Phase diagram for U-Pu-Zr  (Fphase.for — Karahan 2007)
//
// Returns:
//   phase  — "alpha", "beta", "gamma", "delta",
//             "alpha_delta", "beta_gamma"
//   pfrac  — volume fraction of first phase in biphase nodes (0..1)
//
// Inputs:
//   T_K           : fuel temperature [K]
//   pu_weight_frac: Pu weight fraction [-]
//   zr_weight_frac: Zr weight fraction [-]
// -------------------------------------------------------------------------
struct PhaseResult {
    std::string phase;
    double pfrac; // phase fraction (1.0 for single-phase, <1.0 for biphase)
};

inline PhaseResult upuzrPhase(double T_K, double pu_weight_frac, double zr_weight_frac) {
    // Convert weight fractions to molar fractions (same as Fphase.for)
    const double pu = pu_weight_frac;
    const double zr = zr_weight_frac;
    const double ur = 1.0 - pu - zr;

    const double ma = 1.0 / (ur / 238.02891 + pu / 244.06 + zr / 91.22);
    const double zrmolar = (zr / 91.22) * ma;
    const double pumolar  = (pu / 244.06) * ma;
    (void)pumolar; // used in interpolation functions via pu_weight_frac

    // Phase diagram boundary lines (interpolated between pure U-Zr and 19% Pu-Zr)
    // line1: alpha/alpha_delta boundary (lower Zr limit of alpha_delta)
    const double line1_zr     = 0.01;
    const double line1_19pu   = 0.001 + (T_K - 773.15) / 2968.8;
    const double line1 = (1.0 - pu * 100.0 / 19.0) * line1_zr
                       + (pu * 100.0 / 19.0) * line1_19pu;

    // line2: alpha_delta/delta boundary
    const double line2_zr   = (T_K - 813.15) / (935.15 - 813.15) * (0.588 - 0.676) + 0.676;
    const double line2_19pu = 0.539 - (T_K - 773.15) / 9500.0;
    const double line2 = (1.0 - pu * 100.0 / 19.0) * line2_zr
                       + (pu * 100.0 / 19.0) * line2_19pu;

    // line3: alpha/beta boundary temperature
    const double line3_zr   = 935.15;
    const double line3_19pu = 868.15;
    const double line3 = (1.0 - pu * 100.0 / 19.0) * line3_zr
                       + (pu * 100.0 / 19.0) * line3_19pu;

    // line4: beta/beta_gamma lower boundary
    const double line4_zr   = 0.01;
    const double line4_19pu = 0.032 - (T_K - 868.15) / 6111.1;
    const double line4 = (1.0 - pu * 100.0 / 19.0) * line4_zr
                       + (pu * 100.0 / 19.0) * line4_19pu;

    // line5: beta_gamma/gamma boundary
    const double line5_zr   = (T_K - 935.15) / (965.15 - 935.15) * (0.444 - 0.588) + 0.588;
    double line5_19pu;
    if (T_K < 905.0)
        line5_19pu = 0.529 - (T_K - 868.15) / 440.5;
    else
        line5_19pu = 0.445 - (T_K - 905.15) / 200.0;
    const double line5 = (1.0 - pu * 100.0 / 19.0) * line5_zr
                       + (pu * 100.0 / 19.0) * line5_19pu;

    // line6: beta/gamma transition temperature
    const double line6_zr   = 965.15;
    const double line6_19pu = 923.15;
    const double line6 = (1.0 - pu * 100.0 / 19.0) * line6_zr
                       + (pu * 100.0 / 19.0) * line6_19pu;

    // Phase conditions
    const bool cond_alpha       = (zrmolar < line1) && (T_K < line3);
    const bool cond_beta        = (T_K > line3) && (zrmolar < line4) && (T_K < line6);
    const bool cond_alpha_delta = (zrmolar > line1) && (T_K < line3) && (zrmolar < line2);
    const bool cond_beta_gamma  = (T_K > line3) && (T_K < line6)
                                && (zrmolar > line4) && (zrmolar < line5);
    const bool cond_gamma       = (T_K > line6)
                                || ((T_K < line6) && (zrmolar > line5) && (T_K > line3));
    const bool cond_delta       = (T_K < line3) && (zrmolar > line2);

    PhaseResult result;
    result.pfrac = 1.0; // single phase default

    if (cond_gamma)       result.phase = "gamma";
    else if (cond_delta)  result.phase = "delta";
    else if (cond_beta_gamma) {
        result.phase = "beta_gamma";
        result.pfrac = 0.5;
    } else if (cond_beta)  result.phase = "beta";
    else if (cond_alpha_delta) {
        result.phase = "alpha_delta";
        result.pfrac = 0.5;
    } else if (cond_alpha) result.phase = "alpha";
    else                   result.phase = "gamma"; // fallback

    return result;
}

// -------------------------------------------------------------------------
// Sodium infiltration fraction  (Sinf.for)
//
// Returns psod: fraction of pores filled with sodium [0..1].
//
//   bup_FIMA    : local burnup [FIMA]
//   buhard_FIMA : burnup at moment of hard contact [FIMA]
//   radnode     : current radius of fuel node [m]
//   maxradnode  : outer fuel radius [m]
//   phase       : fuel phase string ("alpha", "alpha_delta", etc.)
//   gap_flag    : gap state ("open", "soft", "clos")
// -------------------------------------------------------------------------
inline double sodiumInfiltration(double bup_FIMA, double buhard_FIMA,
                                  double radnode, double maxradnode,
                                  const std::string& phase,
                                  const std::string& gap_flag)
{
    // Only applicable in alpha or alpha_delta phase
    if (phase != "alpha" && phase != "alpha_delta") return 0.0;
    // Only in outer 40% of fuel radius
    if (radnode < 0.6 * maxradnode) return 0.0;

    if (gap_flag == "soft" || gap_flag == "open") {
        return 0.6;
    } else { // "clos"
        return std::max(0.3, 0.6 - 5.0 * (bup_FIMA - buhard_FIMA));
    }
}

// -------------------------------------------------------------------------
// Fuel axial anisotropy factor (fanis.for)
//
// Metallic (U-Pu-Zr) fuel swells anisotropically along its axial grain
// structure. This produces localized fuel-clad contact ("soft" contact,
// see Baseir.for) well before the bulk/average gap has geometrically
// closed. fanis_coef gives the fraction of the as-fabricated gap consumed
// by this local protrusion.
//
//   F    : peak linear-heat-rate parameter [W/cm^2], = ql/(2*rfo0)*1e-4
//   cont : fuel Pu weight fraction [-]
// -------------------------------------------------------------------------
inline double upuzrFanis(double F, double cont) {
    double fanis;
    if (F < 700.0) {
        if      (cont < 0.08)  fanis = 0.15*cont/0.08 + 0.45;
        else if (cont >= 0.19) fanis = 0.62;
        else                   fanis = 0.02*cont/0.11 + 0.60;
    } else if (F >= 900.0) {
        if      (cont < 0.08)  fanis = 0.15*cont/0.08 + 0.45;
        else if (cont >= 0.19) fanis = 0.90;
        else                   fanis = 0.3*cont/0.11 + 0.60;
    } else {
        if      (cont < 0.08)  fanis = 0.15*cont/0.08 + 0.45;
        else if (cont >= 0.19) fanis = 0.62 + 0.28*(F-700.0)/200.0;
        else                   fanis = cont/0.11*(0.02 + 0.28*(F-700.0)/200.0) + 0.6;
    }
    return fanis;
}

// -------------------------------------------------------------------------
// Fission gas release fraction for U-Pu-Zr (empirical Karahan 2009 — FGR.for)
//
//   bupave_FIMA    : average fuel burnup [FIMA] (for threshold check)
//   bup_FIMA_local : local burnup [FIMA]
//
// Returns fgr fraction in [0..1].
// -------------------------------------------------------------------------
inline double upuzrFGRFraction(double bupave_FIMA, double bup_FIMA_local) {
    const double buperave = bupave_FIMA * 100.0; // convert to at%
    if (buperave < 0.8) return 0.0;
    const double buper = bup_FIMA_local * 100.0;
    return std::min(0.8 * (1.0 - std::exp(-buper / 1.8)), 1.0);
}

// -------------------------------------------------------------------------
// Burnup conversion: MWd/kgU → FIMA  (Baseir.for formula, Daniele variant)
//
//   qqv    : power density [W/m3]
//   dt     : time step [s]
//   rof0   : initial fuel density [kg/m3]
//   pu     : Pu weight fraction [-]
//   zr     : Zr weight fraction [-]
// -------------------------------------------------------------------------
inline double upuzrBupFIMAIncrement(double qqv, double dt,
                                     double rof0, double pu, double zr) {
    // Molecular weight of heavy metals only (no Zr counted in HM for FIMA)
    const double molwe = 1.0 / (pu / 0.244 + (1.0 - zr - pu) / 0.238029);
    const double hm = rof0 * MNA_AVOGADRO / molwe; // HM atoms per m3
    // 200 MeV/fiss = 3.204354e-11 J/fiss
    return (qqv * dt) / (3.204354e-11 * hm);
}

// -------------------------------------------------------------------------
// Volumetric swelling strain increment for U-Pu-Zr (simplified solid FP model)
// Based on Ogata (1999) solid FP: 1.5% per at% burnup
//
//   bup_FIMA     : current burnup [FIMA]
//   bup0_FIMA    : burnup at previous step [FIMA]
// Returns volumetric swelling strain increment [-] (to be added to efs).
// -------------------------------------------------------------------------
inline double upuzrSolidSwellingIncrement(double bup_FIMA, double bup0_FIMA) {
    const double dbup_atpct = (bup_FIMA - bup0_FIMA) * 100.0;
    return 0.015 * dbup_atpct; // 1.5% volume per at% (Ogata 1999)
}

// -------------------------------------------------------------------------
// Zirconium redistribution — compute per-node fluxes (Zrdist.for)
// This is a simplified one-step integration for use in the solver post-step.
//
// Takes per-node arrays and updates zr_at, zr_wf, pu_wf, ur_wf in-place.
//
//   nf       : number of fuel radial nodes
//   T_node[] : centre-of-node temperatures [K], size nf (avg of boundaries)
//   rad0[]   : initial node boundary radii [m], size nf+1
//   zr_at[], pu_at[], ur_at[] : per-node atomic fractions (updated in-place)
//   zr_wf[], pu_wf[], ur_wf[] : per-node weight fractions (updated in-place)
//   mass[]   : node mass [kg] (fixed)
//   dvol[]   : node volume [m3] (fixed)
//   c_zr[]   : per-node Zr atomic density [atom/m3] (updated in-place)
//   phase[]  : per-node phase strings (size nf-1, for centre nodes)
//   pfrac[]  : per-node phase fraction (size nf-1)
//   zr_cont  : nominal Zr weight fraction (global fuel composition)
//   pu_cont  : nominal Pu weight fraction
//   dt       : time step [s]
// -------------------------------------------------------------------------
inline void upuzrZirconiumRedistribution(
    int nf,
    const double* T_node,    // size nf-1 (centre temperatures)
    const double* rad0,      // size nf+1 (boundary radii)
    const std::string* phase_node, // size nf-1
    const double* pfrac_node,      // size nf-1
    double* zr_at, double* pu_at, double* ur_at,
    double* zr_wf, double* pu_wf, double* ur_wf,
    double* zr_atoms, double* pu_atoms, double* ur_atoms,
    double* c_zr, const double* mass, const double* dvol,
    double dt)
{
    const int nc = nf - 1; // number of centre nodes
    const double R_kJ = 8.314e-3; // kJ/(mol·K)
    const double lim_alpha = 0.05, lim_delta = 0.40, lim_beta = 0.05;

    // Compute per-node diffusion coefficients and heat of transport
    std::vector<double> Dc1(nc), Dc2(nc), Qc1(nc), Qc2(nc);
    for (int i = 0; i < nc; ++i) {
        const double tc = T_node[i];
        const double pu_at_i = pu_at[i];
        double D0a, D0d, D0b, D0g, Q0a, Q0d, Q0b, Q0g, Qa_h, Qd_h, Qb_h, Qg_h;
        // See Table 3.1 Timpano thesis: https://www.research-collection.ethz.ch/bitstreams/7f277774-0833-4889-9bda-82a068533c13/download
        if (pu_at_i < 0.07) {
            D0a = 2e-7; D0d = 2e-7; D0b = 5.7e-6;
            D0g = 0.4 * std::pow(10.0, -5.1 - 8.05*zr_at[i] + 9.13*zr_at[i]*zr_at[i]);
            Q0a = 170; Q0d = 150; Q0b = 180;
            Q0g = 128 - 107*zr_at[i] + 174*zr_at[i]*zr_at[i];
            Qa_h = 0; Qd_h = 0; Qb_h = 0; Qg_h = -80;
        } else {
            D0a = 4e-6; D0d = 4e-6; D0b = 7e-4;
            D0g = std::pow(10.0, -5.1 - 8.05*zr_at[i] + 9.13*zr_at[i]*zr_at[i]);
            Q0a = 170; Q0d = 150; Q0b = 180;
            Q0g = 128 - 107*zr_at[i] + 174*zr_at[i]*zr_at[i];
            Qa_h = 150; Qd_h = 150; Qb_h = 650; Qg_h = -200;
        }
        const double Dea = D0a * std::exp(-Q0a / (R_kJ * tc));
        const double Ded = D0d * std::exp(-Q0d / (R_kJ * tc));
        const double Deb = D0b * std::exp(-Q0b / (R_kJ * tc));
        const double Deg = D0g * std::exp(-Q0g / (R_kJ * tc));
        const auto& ph = phase_node[i];
        if (ph == "alpha")       { Dc1[i]=Dea; Dc2[i]=0; Qc1[i]=Qa_h; Qc2[i]=0; }
        else if (ph == "delta")  { Dc1[i]=Ded; Dc2[i]=0; Qc1[i]=Qd_h; Qc2[i]=0; }
        else if (ph == "beta")   { Dc1[i]=Deb; Dc2[i]=0; Qc1[i]=Qb_h; Qc2[i]=0; }
        else if (ph == "gamma")  { Dc1[i]=Deg; Dc2[i]=0; Qc1[i]=Qg_h; Qc2[i]=0; }
        else if (ph == "alpha_delta") { Dc1[i]=Dea; Dc2[i]=Ded; Qc1[i]=Qa_h; Qc2[i]=Qd_h; }
        else if (ph == "beta_gamma")  { Dc1[i]=Deb; Dc2[i]=Deg; Qc1[i]=Qb_h; Qc2[i]=Qg_h; }
        else                     { Dc1[i]=Dea; Dc2[i]=0; Qc1[i]=0; Qc2[i]=0; }
    }

    // Subdivide timestep for convergence (20 sub-steps as in Fortran)
    const int ntime = 20;
    const double dtnew = dt / ntime;

    std::vector<double> c_zr_old(nc);
    for (int i = 0; i < nc; ++i) c_zr_old[i] = c_zr[i];

    std::vector<double> Jdi(nc+1, 0.0), Jtr(nc+1, 0.0);

    // Temperature at each centre node (T_node) for T gradient
    // rad boundaries: rad0[0..nc] (nf-1 boundaries between nf nodes)

    for (int t = 0; t < ntime; ++t) {
        Jdi[0] = 0; Jtr[0] = 0; Jdi[nc] = 0; Jtr[nc] = 0;
        for (int i = 1; i < nc; ++i) {
            const double dT = T_node[i] - T_node[i-1];
            const double dr = rad0[i] - rad0[i-1];
            const double c_mid = 0.5 * (c_zr_old[i] + c_zr_old[i-1]);
            const double Dm1   = 0.5 * (Dc1[i] + Dc1[i-1]); // Zr diffusion coefficient [m2/s]
            const double Qm1   = 0.5 * (Qc1[i] + Qc1[i-1]); // Q heat of transport of local phase [J/mol]
            const double Tm    = 0.5 * (T_node[i] + T_node[i-1]);

            const auto& ph = phase_node[i];
            if (ph == "alpha_delta" || ph == "beta_gamma") {
                // Biphase: no diffusive current, sum transport currents from each phase
                Jdi[i] = 0;
                const double pf = pfrac_node[i];
                const double Dm2 = 0.5 * (Dc2[i] + Dc2[i-1]);
                const double Qm2 = 0.5 * (Qc2[i] + Qc2[i-1]);
                Jtr[i] = -pf * dT / dr * (Dm1 * Qm1) / (R_kJ * Tm * Tm) * c_mid
                        - (1.0 - pf) * dT / dr * (Dm2 * Qm2) / (R_kJ * Tm * Tm) * c_mid;
            } else {
                Jdi[i] = -Dm1 * (c_zr_old[i] - c_zr_old[i-1]) / dr;
                Jtr[i] = -dT / dr * (Dm1 * Qm1) / (R_kJ * Tm * Tm) * c_mid;
            }
        }

        // Cut currents at solubility limits
        for (int i = 0; i < nc; ++i) {
            const auto& ph = phase_node[i];
            if (zr_at[i] < lim_alpha && (ph == "alpha" || ph == "alpha_delta")) {
                Jtr[i+1] = 0; Jdi[i+1] = 0;
            } else if (zr_at[i] < lim_beta &&
                       (ph == "beta" || ph == "beta_gamma" || ph == "gamma")) {
                Jtr[i] = 0; Jtr[i+1] = 0; Jdi[i+1] = 0;
            } else if (zr_at[i] >= lim_delta &&
                       (ph == "beta_gamma" || ph == "gamma" || ph == "alpha_delta")) {
                Jtr[i] = 0; Jtr[i+1] = 0; Jdi[i+1] = 0;
            }
        }

        // Update Zr atomic density
        for (int i = 0; i < nc; ++i) {
            const double r0 = rad0[i], r1 = rad0[i+1];
            c_zr[i] = c_zr_old[i] + 2.0 * dtnew / (r1*r1 - r0*r0)
                     * (r0 * (Jdi[i] + Jtr[i]) - r1 * (Jdi[i+1] + Jtr[i+1]));
        }

        // Recompute atomic fractions
        for (int i = 0; i < nc; ++i) {
            zr_atoms[i] = c_zr[i] * dvol[i];
            ur_atoms[i] = (mass[i] * MNA_AVOGADRO * 1e3
                          - (zr_atoms[i] * 91.22 + pu_atoms[i] * 244.06)) / 238.02891;
            const double tot = zr_atoms[i] + ur_atoms[i] + pu_atoms[i];
            zr_at[i] = zr_atoms[i] / tot;
            ur_at[i] = ur_atoms[i] / tot;
            pu_at[i] = pu_atoms[i] / tot;
            const double ma_v = zr_at[i]*91.22 + ur_at[i]*238.02891 + pu_at[i]*244.06;
            zr_wf[i] = zr_at[i] * 91.22 / ma_v;
            pu_wf[i] = pu_at[i] * 244.06 / ma_v;
            ur_wf[i] = ur_at[i] * 238.02891 / ma_v;
            c_zr_old[i] = c_zr[i];
        }
    }
}

// -------------------------------------------------------------------------
// Cladding wastage (Clanth.for — lanthanide diffusion model for HT-9)
//
//   nz       : number of axial layers
//   dt       : time step [s]
//   time_s   : current time [s]  (> 0)
//   tstep    : time step index (1-based; initialise on first call)
//   qqv[]    : power density per layer [W/m3], size nz
//   cit[]    : cladding inner temperature per layer [K], size nz
//   flag[]   : gap state strings per layer, size nz
//   rof0     : initial fuel density [kg/m3]
//   pucont, zrcont : fuel composition
//   rfo0, dz0      : fuel geometry
//   ntot[]   : fuel atom density per layer [atom/m3] (initialised on tstep==1)
//   clfuel[] : lanthanide concentration per layer [-] (updated)
//   xwast[]  : wastage thickness per layer [m] (updated)
// -------------------------------------------------------------------------
inline void ht9CladWastage(
    int nz, double dt, double time_s,
    const double* qqv, const double* cit, const std::string* flag,
    double rof0, double pucont, double zrcont, double rfo0, const double* dz0,
    double* ntot, double* clfuel, double* xwast)
{
    // Model parameters from Clanth.for (fitted to SAS data)
    const double csolfuel = 0.0;     // lanthanide solubility in fuel [-]
    const double csolclad = 0.1;     // lanthanide solubility in cladding [-]
    const double qlanth   = 333962.4;// activation energy [J/mol]
    const double d0l      = 1632573.0;// diffusion pre-exp [m2/s]
    const double yl       = 15.29e9; // fission yield of lanthanides [1/J]
    const double R_Jmol   = 8.314;   // J/(mol·K)

    // Average atomic mass of fuel [kg/mol]
    const double ur = 1.0 - pucont - zrcont;
    const double ma = 1.0 / (ur / 238.02891e-3 + pucont / 244.06e-3 + zrcont / 91.22e-3);

    const int ntime = 20;
    const double dtnew = dt / ntime;

    for (int t = 0; t < ntime; ++t) {
        for (int j = 0; j < nz; ++j) {
            // Self-initialise on first call (ntot==0 means uninitialised).
            if (ntot[j] <= 0.0) {
                ntot[j]   = rof0 * MNA_AVOGADRO / ma;
                clfuel[j] = 0.0;
                xwast[j]  = 0.0;
            }
            // Diffusion coefficient
            double dl;
            const double Tcit = std::max(cit[j], 800.0);
            dl = d0l * std::exp(-qlanth / (R_Jmol * Tcit));

            // Update lanthanide concentration in fuel
            clfuel[j] += (yl * qqv[j] * dtnew) / ntot[j];

            // Wastage only during soft or closed contact
            if (flag[j] == "soft" || flag[j] == "clos") {
                xwast[j] += 0.5 * (clfuel[j] - csolfuel) / (csolclad - csolfuel)
                           * std::sqrt(dl / std::max(time_s, 1.0)) * dtnew;
            } else {
                xwast[j] = 0.0;
            }
        }
    }
}

} // namespace fred
