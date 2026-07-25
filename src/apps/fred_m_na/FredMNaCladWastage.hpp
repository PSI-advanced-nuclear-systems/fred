#pragma once
// Cladding wastage models for U-Pu-Zr fuel in HT-9 cladding.
//
// Two selectable models (FredMNaSolver::setCladWastageModel):
//
//   PrecipitationKinetics (default) — layer-lumped lanthanide inventory with
//     the sqrt(D/t) precipitation-front law.  Direct port of Clanth.for
//     (Timpano/EPFL FRED_M_OCT24.SRC), coefficients refitted to SAS4A data.
//
//   LaTracking — per-fuel-node lanthanide diffusion following SAS4A/MFUEL
//     (SAS4A/SASSYS-1 v5.7 Theory Manual, Appendix 9.3; see
//     issues/lanthanide_diffusion_mfuel.md):
//       du/dt = (1/r) d/dr ( r D_La(T) du/dr ) + Y_La q'''
//     with symmetry at r=0 and, while the gap is in soft/clos contact, a
//     Dirichlet sink u(R_fc)=0 at the fuel-cladding interface (same contact
//     gating as the PK model; with an open gap the outer face is zero-flux
//     and no wastage grows).  The wastage depth integrates the arriving
//     interface flux:  delta = (1/u_sol) * integral J_La dt.
//
// Neither model resets accumulated wastage: the open-gap reset branch of
// legacy Clanth.for is unreachable for U-Pu-Zr (the gap-contact ratchet
// never reopens; Baseir.for's reopen check is gated fmat.ne.'upuzr').

#include "apps/fred_m_na/FredMNaCladdingMaterial.hpp"
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>

namespace fred {

enum class CladWastageModel {
    PrecipitationKinetics = 1,   // Clanth.for port (SAS-refitted coefficients)
    LaTracking            = 2    // MFUEL App. 9.3 radial diffusion model
};

// ---------------------------------------------------------------------------
// Model parameters — simulation/scheme-level knobs only.  Per-material
// correlation values (diffusion coefficients, solubilities, uptake
// capacity) live on the FredMNaCladdingMaterial passed to the functions
// below (wastageDiffusionPK, wastageSolubilityFuel/Clad,
// wastageDiffusionLaTracking, wastageUptakeCapacity), not here.
// ---------------------------------------------------------------------------
struct CladWastageParams {
    // Shared
    double y_la     = 15.29e9;    // lanthanide fission yield [atoms/J]
    double T_floor  = 800.0;      // irradiation-enhanced diffusivity floor [K]

    // true  (default): interface sink only during soft/clos contact
    //                  (same gating as the PK model)
    // false: sink active from t=0, as in MFUEL's unconditional Dirichlet BC
    //        (lanthanide transport through the sodium bond) — gives the
    //        flatter axial wastage profile MFUEL shows at cold layers
    bool la_sink_requires_contact = true;
    // Radial subgrid refinement for the La diffusion solve.  The floored
    // diffusivity makes the interface boundary layer ~sqrt(2 D t) ~ 0.3 mm,
    // thinner than a typical thermal cell (0.39 mm at nf=11), so solving on
    // the thermal grid under-resolves the sink flux (mesh study: 28 -> 48 um
    // at 1825 d going from nf=11 to an 81-point subgrid, IF-AVG layer 12).
    // The solve runs on (nf-1)*la_refine + 1 uniform subnodes with linearly
    // interpolated temperatures.
    int la_refine = 8;
};

// ---------------------------------------------------------------------------
// PrecipitationKinetics (Clanth.for — lanthanide precipitation model, HT-9)
//
//   nz       : number of axial layers
//   dt       : time step [s]
//   time_s   : current time [s]  (> 0)
//   qqv[]    : power density per layer [W/m3], size nz
//   cit[]    : cladding inner temperature per layer [K], size nz
//   flag[]   : gap state strings per layer, size nz
//   rof0     : initial fuel density [kg/m3]
//   pucont, zrcont : fuel composition
//   rfo0, dz0      : fuel geometry
//   ntot[]   : fuel atom density per layer [atom/m3] (initialised on first call)
//   clfuel[] : lanthanide concentration per layer [-] (updated)
//   xwast[]  : wastage thickness per layer [m] (updated)
//   clad     : cladding material — supplies wastageDiffusionPK(T),
//              wastageSolubilityFuel(), wastageSolubilityClad()
// -------------------------------------------------------------------------
inline void cladWastagePrecipitation(
    int nz, double dt, double time_s,
    const double* qqv, const double* cit, const std::string* flag,
    double rof0, double pucont, double zrcont, double rfo0, const double* dz0,
    double* ntot, double* clfuel, double* xwast,
    const FredMNaCladdingMaterial& clad,
    const CladWastageParams& p = CladWastageParams{})
{
    (void)rfo0; (void)dz0;
    constexpr double AVO    = 6.02214076e23;

    // Average atomic mass of fuel [kg/mol]
    const double ur = 1.0 - pucont - zrcont;
    const double ma = 1.0 / (ur / 238.02891e-3 + pucont / 244.06e-3 + zrcont / 91.22e-3);

    const int ntime = 20;
    const double dtnew = dt / ntime;
    const double csol_fuel = clad.wastageSolubilityFuel();
    const double csol_clad = clad.wastageSolubilityClad();

    for (int t = 0; t < ntime; ++t) {
        for (int j = 0; j < nz; ++j) {
            // Self-initialise on first call (ntot==0 means uninitialised).
            if (ntot[j] <= 0.0) {
                ntot[j]   = rof0 * AVO / ma;
                clfuel[j] = 0.0;
                xwast[j]  = 0.0;
            }
            const double Tcit = std::max(cit[j], p.T_floor);
            const double dl = clad.wastageDiffusionPK(Tcit);

            // Update lanthanide concentration in fuel
            clfuel[j] += (p.y_la * qqv[j] * dtnew) / ntot[j];

            // Wastage only during soft or closed contact
            if (flag[j] == "soft" || flag[j] == "clos") {
                xwast[j] += 0.5 * (clfuel[j] - csol_fuel)
                           / (csol_clad - csol_fuel)
                           * std::sqrt(dl / std::max(time_s, 1.0)) * dtnew;
            } else {
                // Unreachable after first contact for upuzr (monotonic
                // gap ratchet); zeroes only the pre-contact value (= 0).
                xwast[j] = 0.0;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// LaTracking — per-node radial lanthanide diffusion (MFUEL App. 9.3).
//
// One axial layer per call.  Explicit conservative FVM on the vertex-centred
// fuel node grid (uniform spacing dr = rfo0/(nf-1), same grid as
// HeatConduction), internally sub-stepped to the diffusion stability limit.
//
//   nf      : number of radial fuel nodes
//   dt      : time step [s]
//   qqv     : layer power density [W/m3] (source is radially uniform)
//   T_fuel  : fuel node temperatures [K], size nf
//   T_ci    : cladding inner temperature [K] — evaluates the surface-node
//             diffusivity: the "final" diffusivity a lanthanide sees before
//             entering the cladding is controlled by the clad temperature
//             (coolant-side cooling suppresses La transport into the clad)
//   rfo0    : as-fabricated fuel outer radius [m]
//   contact : true while gap state is "soft" or "clos"
//   clad    : cladding material — supplies wastageDiffusionLaTracking(T),
//             wastageUptakeCapacity()
//   p       : simulation/scheme parameters (subgrid refinement, T floor, ...)
//   c_sub   : subgrid La number density [#/m3], self-sizing to
//             (nf-1)*p.la_refine + 1 (persistent state, checkpointed)
//   c_la    : OUT — La density sampled back at the nf fuel nodes (for H5)
//   xwast   : wastage thickness [m] (updated; monotonic)
//
// Wastage bookkeeping is done by exact inventory conservation: per sub-step,
// atoms absorbed by the sink = generated - d(inventory), which equals the
// time-integrated interface flux without needing a one-sided gradient.
// ---------------------------------------------------------------------------
inline void cladWastageLaTracking(
    int nf, double dt, double qqv,
    const double* T_fuel, double T_ci, double rfo0, bool contact,
    const FredMNaCladdingMaterial& clad,
    const CladWastageParams& p,
    std::vector<double>& c_sub, double* c_la, double& xwast)
{
    const int m  = std::max(1, p.la_refine);
    const int n  = (nf - 1) * m + 1;         // subgrid node count
    if ((int)c_sub.size() != n) c_sub.assign(n, 0.0);
    const double dr = rfo0 / (n - 1);
    const double S  = p.y_la * qqv;          // source [#/m3/s]

    // Subgrid temperatures: linear interpolation between fuel nodes; the
    // surface node uses the cladding inner temperature (rate-limiting
    // "final" diffusivity).  Diffusion coefficients per subnode.
    std::vector<double> D(n);
    double Dmax = 0.0;
    for (int k = 0; k < n; ++k) {
        double T;
        if (k == n - 1) {
            T = T_ci;
        } else {
            const int    i = k / m;
            const double w = double(k - i * m) / m;
            T = (1.0 - w) * T_fuel[i] + w * T_fuel[std::min(i + 1, nf - 1)];
        }
        T = std::max(T, p.T_floor);
        D[k] = clad.wastageDiffusionLaTracking(T);
        Dmax = std::max(Dmax, D[k]);
    }

    // Sub-step to the explicit stability limit (factor 0.25 safety; the
    // centre-node update has the stiffest coefficient 4D/dr^2).
    const double dt_stab = 0.25 * dr * dr / std::max(4.0 * Dmax, 1e-30);
    const int n_sub = std::max(1, (int)std::ceil(dt / dt_stab));
    const double h = dt / n_sub;

    // Control volumes (per unit length / 2*pi): vertex-centred cells.
    //   V_0 = (dr/2)^2/2, V_k = r_k*dr (interior), V_last = half cell at edge
    // Faces at r_{k+1/2} = (k+0.5)*dr with arithmetic-mean D.
    std::vector<double> V(n);
    V[0] = 0.5 * (0.5 * dr) * (0.5 * dr);
    for (int k = 1; k < n - 1; ++k) V[k] = (k * dr) * dr;
    V[n-1] = ((n - 1) * dr - 0.25 * dr) * (0.5 * dr);

    // SIGN CONVENTION: OUTWARD IS POSITIVE.  flux[k] is the (r-weighted)
    // diffusive transport through face k+1/2 in the +r direction:
    //   flux[k] = r_{k+1/2} * D_f * (c[k] - c[k+1]) / dr   [= -r*D*du/dr]
    // positive when atoms flow outward (node k richer than node k+1).
    // Node balances read (flux[k-1] - flux[k])/V: gain through the inner
    // face, loss through the outer face.
    std::vector<double> flux(n - 1);
    for (int s = 0; s < n_sub; ++s) {
        for (int k = 0; k < n - 1; ++k) {
            const double rf  = (k + 0.5) * dr;
            const double Df  = 0.5 * (D[k] + D[k+1]);
            flux[k] = rf * Df * (c_sub[k] - c_sub[k+1]) / dr;
        }

        // Inventory before (for conservation bookkeeping)
        double inv0 = 0.0;
        for (int k = 0; k < n; ++k) inv0 += V[k] * c_sub[k];

        // Centre node: symmetry (no inner face), loses through its outer face
        c_sub[0] += h * (-flux[0] / V[0] + S);
        for (int k = 1; k < n - 1; ++k)
            c_sub[k] += h * ((flux[k-1] - flux[k]) / V[k] + S);

        // Outer node: half control volume at the interface.  Open gap:
        // zero-flux outer face.  Contact: Dirichlet u=0 applied at the
        // boundary FACE (r = rfo0), so the sink drains the surface node
        // through the final face diffusivity D[n-1] = D(T_ci) over dr/2 —
        // the cladding temperature is rate-limiting, and the node's own
        // fission source is not spuriously dumped into the sink.
        double sink_flux = 0.0;   // outward flux into the clad sink
        if (contact) {
            const double rR = (n - 1) * dr;   // = rfo0
            // (c_sub[n-1] - 0.0): gradient toward the boundary-face value
            // u(r_fo) = 0 over the half-spacing dr/2
            sink_flux = rR * D[n-1] * (c_sub[n-1] - 0.0) / (0.5 * dr);
        }
        c_sub[n-1] += h * ((flux[n-2] - sink_flux) / V[n-1] + S);
        for (int k = 0; k < n; ++k) c_sub[k] = std::max(c_sub[k], 0.0);

        if (contact) {
            // absorbed = generated - d(inventory): exactly the sink-face
            // outflow for this conservative scheme (per unit length / 2*pi)
            double inv1 = 0.0, Vtot = 0.0;
            for (int k = 0; k < n; ++k) { inv1 += V[k] * c_sub[k]; Vtot += V[k]; }
            const double absorbed = S * h * Vtot - (inv1 - inv0);
            if (absorbed > 0.0)
                // divide by interface circumference (rfo0, per unit length /
                // 2*pi) and by the clad uptake capacity -> wastage depth
                xwast += absorbed / (rfo0 * clad.wastageUptakeCapacity());
        }
    }

    // Sample the subgrid back onto the nf fuel nodes (state/H5 output)
    for (int i = 0; i < nf; ++i) c_la[i] = c_sub[i * m];
}

} // namespace fred
