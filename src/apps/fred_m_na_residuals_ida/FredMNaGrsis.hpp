#pragma once
// GRSIS bubble-swelling and fission-gas-release model for U-Pu-Zr metallic fuel.
//
// Ported from GRSIS_METALFUEL.for (Timpano/EPFL FRED_M_OCT24.SRC).
// Reference: Karahan (2009), Ch. 4. GRSIS model originally from Lee (1999).
//
// Two data-mode presets are available (FEAST is the default in the Fortran):
//   FEAST — calibrated to FEAST-METAL code predictions
//   GRSIS  — original GRSIS parameters
//
// Usage (call once per fuel node per time step):
//   upuzrGrsis(state, first_step, dt, T_K, fissd, bup_FIMA, bup0_FIMA,
//              gpres_Pa, pfc_Pa, sigfr_Pa, sigfh_Pa, sigfz_Pa,
//              gap_closed, params);
//
// Outputs in GrsisFuelNodeState:
//   swtot   — total swelling fraction [-] (solid FP + gas bubbles)
//   swopen  — open-bubble swelling fraction [-]  (= sw_v3)
//   swclose — closed-bubble swelling fraction [-]
//   swgas   — total gas swelling [-]
//   swsol   — solid-FP swelling [-]
//   frtot   — gas-filled porosity [-]
//   cgrel   — cumulative released gas concentration [atom/m3]
//   cggen   — cumulative generated gas concentration [atom/m3]

#include <cmath>
#include <algorithm>

namespace fred {

// ---------------------------------------------------------------------------
// Data-mode presets (correspond to the 'dataMode' string in the Fortran)
// ---------------------------------------------------------------------------
enum class GrsisDataMode { FEAST, GRSIS };

struct GrsisParams {
    double rbb1;    // initial radius of small (type-1) bubbles [m]
    double rbb2;    // initial radius of large (type-2) bubbles [m]
    double surt;    // surface tension [N/m]
    double dgo;     // gas diffusion pre-exponential [m2/s]
    double qg;      // gas diffusion activation energy [cal/mol]
    double dso;     // surface diffusion pre-exponential [m2/s]
    double qss;     // surface diffusion activation energy [cal/mol]
    double aa;      // area per surface molecule [m2]
    double kbnucl;  // bubble-1 nucleation constant [bub/(s·atom)]
    double egb1;    // bias factor, gas diffusion to closed bubble-1
    double egb2;    // bias factor, gas diffusion to closed bubble-2
    double egb3;    // bias factor, gas diffusion to open bubble-3
    double ebbc;    // bias factor, bubble diffusion to closed bubble
    double ebbo;    // bias factor, bubble diffusion to open bubble
    double sth;     // threshold closed-bubble swelling for interconnection [-]
    double fth;     // fraction of bubbles interconnected at threshold [-]
    double fv;      // correction for bubble volume after becoming open
    double fs;      // correction for bubble area after becoming open
    double dd1;     // open-bubble formation coefficient [-]

    static GrsisParams feast() {
        return {
            0.5e-6,   // rbb1
            10.0e-6,  // rbb2
            0.8,      // surt
            2.3e-3,   // dgo
            52000.0,  // qg
            2.3,      // dso
            52000.0,  // qss
            9.0e-20,  // aa
            1.0e-20,  // kbnucl
            1.0, 1.0, 1.0, // egb1/2/3
            1.0, 1.0,      // ebbc, ebbo
            0.1,      // sth
            0.01,     // fth
            1.0,      // fv
            0.2,      // fs
            0.235     // dd1
        };
    }

    static GrsisParams grsis() {
        return {
            0.5e-6,    // rbb1
            12.5e-6,   // rbb2
            1.0,       // surt
            0.95e-8,   // dgo
            32000.0,   // qg
            0.95e-5,   // dso
            32000.0,   // qss
            9.0e-20,   // aa
            1.0e-20,   // kbnucl
            1.0, 1.0, 1.0,
            1.0, 1.0,
            0.2,       // sth
            0.3,       // fth
            1.0,       // fv
            0.6,       // fs
            0.235      // dd1
        };
    }
};

// ---------------------------------------------------------------------------
// Per-fuel-node GRSIS bubble state
// ---------------------------------------------------------------------------
struct GrsisFuelNodeState {
    // Gas concentrations [atom/m3]
    double cgfm  = 0.0;  // gas atoms in fuel matrix
    double cgb1  = 0.0;  // gas in closed small bubbles
    double cgb2  = 0.0;  // gas in closed large bubbles
    double cgb31 = 0.0;  // gas in open small bubbles
    double cgb32 = 0.0;  // gas in open large bubbles
    double cggen = 0.0;  // cumulative generated [atom/m3]
    double cgrel = 0.0;  // cumulative released [atom/m3]

    // Bubble dimensions
    double rb1 = 0.0;   // radius of small bubbles [m]
    double rb2 = 0.0;   // radius of large bubbles [m]
    double vb1 = 0.0;   // volume of small bubble [m3]
    double vb2 = 0.0;   // volume of large bubble [m3]

    // Swelling fractions [-]
    double swclose = 0.0;  // by closed bubbles
    double swopen  = 0.0;  // by open bubbles (= sw_v3)
    double swgas   = 0.0;  // total gas swelling
    double swsol   = 0.0;  // solid fission-product swelling
    double swtot   = 0.0;  // total swelling (must not decrease)

    // Porosity fractions [-]
    double frg   = 0.0;  // gas-filled closed porosity (= swclose)
    double frtot = 0.0;  // total gas porosity (= swgas)
};

// ---------------------------------------------------------------------------
// upuzrGrsis — advance one node by dt seconds (explicit forward Euler)
//
//   s          : current node state (updated in place)
//   first_step : true on first call (initialises bubble radii/volumes)
//   dt         : time step [s]
//   T_K        : node temperature [K]
//   fissd      : fission density rate [fiss/m3/s]
//   bup_FIMA   : burnup at end of step [FIMA]
//   bup0_FIMA  : burnup at start of step [FIMA]
//   gpres_Pa   : gap gas pressure [Pa]
//   pfc_Pa     : pellet-clad contact pressure [Pa]
//   sigfr/h/z  : fuel stress components [Pa] (used only in hot-pressing)
//   gap_closed : true when gap is in 'clos' state (legacy GRSIS_metalfuel.for
//                gates hot-pressing on 'clos' only, never 'soft' — see
//                FredMNaSolver.cpp's `!s.gapOpen` call sites)
//   p          : GRSIS parameter set (use GrsisParams::feast())
// ---------------------------------------------------------------------------
inline void upuzrGrsis(GrsisFuelNodeState& s,
                        bool   first_step,
                        double dt,
                        double T_K,
                        double fissd,
                        double bup_FIMA,
                        double bup0_FIMA,
                        double gpres_Pa,
                        double pfc_Pa,
                        double sigfr_Pa,
                        double sigfh_Pa,
                        double sigfz_Pa,
                        bool   gap_closed,
                        const GrsisParams& p)
{
    using std::max; using std::min; using std::exp; using std::pow; using std::sqrt;
    using std::abs;

    constexpr double KB      = 1.38064852e-23; // Boltzmann constant [J/K]
    constexpr double R_CAL   = 1.987;          // gas constant [cal/(mol·K)]
    constexpr double PI      = 3.14159265358979323846;
    constexpr double TINY    = 1.0e-300;

    // -----------------------------------------------------------------------
    // Initialise on first step
    // -----------------------------------------------------------------------
    if (first_step) {
        s.cgfm  = 0.0;
        s.cgb1  = 0.0; s.cgb2  = 0.0;
        s.cgb31 = 0.0; s.cgb32 = 0.0;
        s.cggen = 0.0; s.cgrel = 0.0;
        s.rb1   = p.rbb1; s.rb2 = p.rbb2;
        s.vb1   = 4.0 / 3.0 * PI * p.rbb1 * p.rbb1 * p.rbb1;
        s.vb2   = 4.0 / 3.0 * PI * p.rbb2 * p.rbb2 * p.rbb2;
        s.swclose = 0.0; s.swopen = 0.0;
        s.swgas  = 0.0; s.swsol  = 0.0; s.swtot = 0.0;
        s.frg    = 0.0; s.frtot  = 0.0;
        return;
    }

    // -----------------------------------------------------------------------
    // Save previous swelling for hot-press and monotonicity constraint
    // -----------------------------------------------------------------------
    const double sw_v3_prev  = s.swopen;   // = sw_v3old in Fortran
    const double swtot_prev  = s.swtot;

    // -----------------------------------------------------------------------
    // Hydrostatic stress: FEAST/GRSIS Fortran convention: shyd = gpres + 2 MPa
    // -----------------------------------------------------------------------
    const double shyd = gpres_Pa + 2.0e6;

    // -----------------------------------------------------------------------
    // Gas density per bubble: rhogi = vbi / (kB*T/(2γ/ri + shyd) + 85e-30)
    // Result: atoms per bubble [atoms/bub].  85e-30 m3 ≈ solid-density floor.
    // -----------------------------------------------------------------------
    const double vol_per_atom1 = KB * T_K / (2.0 * p.surt / max(s.rb1, 1.0e-12) + shyd)
                               + 85.0e-30;
    const double vol_per_atom2 = KB * T_K / (2.0 * p.surt / max(s.rb2, 1.0e-12) + shyd)
                               + 85.0e-30;
    const double rhog1 = max(s.vb1 / max(vol_per_atom1, TINY), 0.0);
    const double rhog2 = max(s.vb2 / max(vol_per_atom2, TINY), 0.0);

    // -----------------------------------------------------------------------
    // Number densities [bub/m3]
    // -----------------------------------------------------------------------
    const double nb1  = (rhog1 > 0.0) ? s.cgb1  / rhog1 : 0.0;
    const double nb2  = (rhog2 > 0.0) ? s.cgb2  / rhog2 : 0.0;
    const double nb31 = (rhog1 > 0.0) ? s.cgb31 / rhog1 : 0.0;
    const double nb32 = (rhog2 > 0.0) ? s.cgb32 / rhog2 : 0.0;

    // -----------------------------------------------------------------------
    // Gas-atom diffusion coefficient [m2/s]  dg = dgo * exp(-qg / (R_cal * T))
    // -----------------------------------------------------------------------
    const double dg = p.dgo * exp(-p.qg / (R_CAL * T_K));

    // -----------------------------------------------------------------------
    // Diffusion rate constants [m3/s]
    // -----------------------------------------------------------------------
    const double kg1 = p.egb1 * 4.0 * PI * s.rb1 * dg;
    const double kg2 = p.egb2 * 4.0 * PI * s.rb2 * dg;

    // -----------------------------------------------------------------------
    // Atomic flux into bubble-i by diffusion, jgi [atom/m3/s]
    // -----------------------------------------------------------------------
    const double jg1  = kg1 * s.cgfm * nb1;
    const double jg2  = kg2 * s.cgfm * nb2;
    const double jg31 = kg1 * s.cgfm * nb31;
    const double jg32 = kg2 * s.cgfm * nb32;

    // -----------------------------------------------------------------------
    // Bubble nucleation rate [atom/m3/s]
    // -----------------------------------------------------------------------
    const double jnucl = p.kbnucl * s.cgfm * rhog1;

    // -----------------------------------------------------------------------
    // Bubble diffusion coefficients [m2/s]
    //   db = 1.5 * aa^2 / (pi * r^4) * 1000 * dg
    // -----------------------------------------------------------------------
    const double r1sq = s.rb1 * s.rb1;
    const double r2sq = s.rb2 * s.rb2;
    const double db1  = 1.5 * p.aa * p.aa / (PI * r1sq * r1sq) * 1000.0 * dg;
    const double db2  = 1.5 * p.aa * p.aa / (PI * r2sq * r2sq) * 1000.0 * dg;

    // -----------------------------------------------------------------------
    // Collision constants [m3/s]
    // -----------------------------------------------------------------------
    const double k11 = p.ebbc * 4.0 * PI * (2.0 * s.rb1) * (2.0 * db1);
    const double k12 = p.ebbc * 4.0 * PI * (s.rb1 + s.rb2) * (db1 + db2);
    const double k22 = p.ebbo * 4.0 * PI * (2.0 * s.rb2) * (2.0 * db2);

    // -----------------------------------------------------------------------
    // Coalescence by bubble diffusion [atom/m3/s]
    // -----------------------------------------------------------------------
    const double ab11  = k11 * nb1 * nb1 * 2.0 * rhog1;
    const double ab12  = k12 * nb1 * nb2 * rhog1;
    const double ab131 = k11 * nb1 * nb31 * rhog1;
    const double ab132 = k12 * nb1 * nb32 * rhog1;
    const double ab231 = k12 * nb2 * nb31 * rhog2;
    const double ab232 = k22 * nb2 * nb32 * rhog2;
    const double f12   = (rhog2 > 0.0) ? 2.0 * rhog1 / rhog2 : 0.0;

    // -----------------------------------------------------------------------
    // Bubble growth rates and radius/volume update
    // -----------------------------------------------------------------------
    const double dotrb1 = (s.cgb1 > 0.0) ? jg1 * s.rb1 / (3.0 * s.cgb1) : 0.0;
    const double dotrb2 = (s.cgb2 > 0.0) ? jg2 * s.rb2 / (3.0 * s.cgb2) : 0.0;
    const double dotvb1 = (s.cgb1 > 0.0) ? jg1 * s.vb1 / s.cgb1         : 0.0;
    const double dotvb2 = (s.cgb2 > 0.0) ? jg2 * s.vb2 / s.cgb2         : 0.0;

    s.rb1 = max(s.rb1 + dotrb1 * dt, p.rbb1);
    s.rb2 = max(s.rb2 + dotrb2 * dt, p.rbb2);
    s.vb1 = max(s.vb1 + dotvb1 * dt, 4.0/3.0*PI*p.rbb1*p.rbb1*p.rbb1);
    s.vb2 = max(s.vb2 + dotvb2 * dt, 4.0/3.0*PI*p.rbb2*p.rbb2*p.rbb2);

    // -----------------------------------------------------------------------
    // Coalescence by radial growth: inter-bubble collision probabilities [1/s]
    // l_i = 1.122 / nb_i^(1/3)  (mean inter-bubble distance)
    // p_ij = dotrb_i / (0.5*l_j - rbb_i - rbb_j)
    // -----------------------------------------------------------------------
    double gab11 = 0.0, gab12 = 0.0, gab21 = 0.0;
    double gab131 = 0.0, gab132 = 0.0;
    double gab231 = 0.0, gab232 = 0.0;

    if (nb1 > 0.0) {
        const double l1   = 1.122 / pow(nb1, 1.0/3.0);
        const double denom11 = 0.5 * l1 - 2.0 * p.rbb1;
        const double denom12_from1 = 0.5 * l1 - p.rbb2 - p.rbb1;
        if (denom11 > 0.0)          gab11 = dotrb1 / denom11 * s.cgb1;
        if (nb2 > 0.0) {
            const double l2 = 1.122 / pow(nb2, 1.0/3.0);
            const double d12 = 0.5 * l2 - p.rbb1 - p.rbb2;
            if (d12 > 0.0) gab12 = dotrb1 / d12 * s.cgb1;
        }
        if (denom12_from1 > 0.0)    gab21 = dotrb2 / denom12_from1 * s.cgb2;
    }
    if (nb31 > 0.0) {
        const double l31  = 1.122 / pow(nb31, 1.0/3.0);
        const double d131 = p.dd1 * l31;
        if (d131 > 0.0) {
            gab131 = dotrb1 / d131 * s.cgb1;
            gab231 = dotrb2 / d131 * s.cgb2;
        }
    }
    if (nb32 > 0.0) {
        const double l32  = 1.122 / pow(nb32, 1.0/3.0);
        const double d132 = p.dd1 * l32;
        if (d132 > 0.0) {
            gab132 = dotrb1 / d132 * s.cgb1;
            gab232 = dotrb2 / d132 * s.cgb2;
        }
    }

    // -----------------------------------------------------------------------
    // Concentration change rates [atom/m3/s]
    // -----------------------------------------------------------------------
    const double dotcggen = 0.25 * fissd;
    const double dotcgfm  = dotcggen - (jg1 + jg2 + jg31 + jg32) - jnucl;

    const double dotcgb1  = jnucl + jg1
                           - ab12 - (ab131 + ab132) - (gab12 + gab21)
                           - (gab131 + gab132) - f12 * (ab11 + gab11);
    const double dotcgb2  = jg2 + ab12 + (gab12 + gab21) + f12 * (ab11 + gab11)
                           - (ab231 + ab232) - (gab231 + gab232);
    const double dotcgb31 = (ab131 + ab132) + (gab131 + gab132);
    const double dotcgb32 = (ab231 + ab232) + (gab231 + gab232);

    // -----------------------------------------------------------------------
    // Threshold conversion: closed → open bubbles when swclose >= sth
    // -----------------------------------------------------------------------
    const double cgconv13 = (s.swclose >= p.sth) ? p.fth * s.cgb1 : 0.0;
    const double cgconv23 = (s.swclose >= p.sth) ? p.fth * s.cgb2 : 0.0;

    // -----------------------------------------------------------------------
    // Apply concentration updates
    // -----------------------------------------------------------------------
    s.cggen += dotcggen * dt;
    s.cgfm  += dotcgfm  * dt;
    s.cgfm   = max(s.cgfm, 0.0);
    s.cgb1   = max(s.cgb1 + dotcgb1 * dt - cgconv13, 0.0);
    s.cgb2   = max(s.cgb2 + dotcgb2 * dt - cgconv23, 0.0);
    s.cgb31 += dotcgb31 * dt + cgconv13;
    s.cgb32 += dotcgb32 * dt + cgconv23;
    s.cgb31  = max(s.cgb31, 0.0);
    s.cgb32  = max(s.cgb32, 0.0);

    // Released gas: gas diffusing into open pores + open bubble growth
    s.cgrel += (jg31 + jg32) * dt
             + dotcgb31 * dt + cgconv13
             + dotcgb32 * dt + cgconv23;
    s.cgrel  = max(s.cgrel, 0.0);

    // -----------------------------------------------------------------------
    // Swelling: closed bubbles
    // -----------------------------------------------------------------------
    const double rhog1_new = (s.vb1 > 0.0) ?
        max(s.vb1 / max(KB*T_K/(2.0*p.surt/max(s.rb1,1e-12)+shyd)+85e-30, TINY), 0.0) : 0.0;
    const double rhog2_new = (s.vb2 > 0.0) ?
        max(s.vb2 / max(KB*T_K/(2.0*p.surt/max(s.rb2,1e-12)+shyd)+85e-30, TINY), 0.0) : 0.0;

    const double v1 = (rhog1_new > 0.0) ? s.cgb1  / rhog1_new * s.vb1 : 0.0;
    const double v2 = (rhog2_new > 0.0) ? s.cgb2  / rhog2_new * s.vb2 : 0.0;

    // open bubble volume (before hot-press override)
    double sw_v3 = 0.0;
    if (rhog1_new > 0.0) sw_v3 += s.cgb31 / rhog1_new * s.vb1;
    if (rhog2_new > 0.0) sw_v3 += s.cgb32 / rhog2_new * s.vb2;

    // -----------------------------------------------------------------------
    // Hot pressing: reduce sw_v3 by creep when gap is closed
    // -----------------------------------------------------------------------
    if (gap_closed && sw_v3 > 0.0) {
        double alphac = 0.0;
        if (sw_v3 >= 0.1) {
            alphac = 1.0 / 6.0;
        } else if (sw_v3 > 0.0) {
            alphac = (1.0 / 6.0) * pow(sw_v3 / 0.1, 1.5);
        }

        const double sout = min(gpres_Pa, pfc_Pa);
        const double seq  = 3.0 * sqrt(3.0 * alphac) * sout;
        const double s_MPa = seq * 1.0e-6;

        double dotopn;
        if (T_K < 923.15) {
            dotopn = (5.0e3 * s_MPa + 6.0 * pow(s_MPa, 4.5)) * exp(-26170.0 / T_K);
        } else {
            dotopn = 0.08 * pow(s_MPa, 3.0) * exp(-14350.0 / T_K);
        }
        dotopn = 3.0 * sqrt(3.0 * alphac) * dotopn * 0.01;

        const double dbup = bup_FIMA - bup0_FIMA;
        if (T_K > 800.0) {
            sw_v3 = sw_v3_prev - min(abs(dotopn * dt), 1.5 * max(dbup, 0.0));
        } else {
            sw_v3 = sw_v3_prev - 1.5 * max(dbup, 0.0);
        }
        sw_v3 = max(sw_v3, 0.0);
    }

    // -----------------------------------------------------------------------
    // Accumulate swelling fractions
    // -----------------------------------------------------------------------
    s.swclose = v1 + v2;
    s.swopen  = sw_v3;
    s.swgas   = s.swclose + s.swopen;
    s.swsol   = 0.015 * bup_FIMA * 100.0;
    s.swtot   = max(s.swsol + s.swgas, swtot_prev);  // monotonically non-decreasing

    s.frg   = s.swclose;
    s.frtot = s.swgas;
}

} // namespace fred
