#pragma once
#include <vector>
#include <string>

namespace fred {

// Runtime state of a single axial layer of the fuel rod.
// All arrays are indexed 0-based and sized at construction.
//
// Differential variables (time derivatives tracked by the DAE solver):
//   T[i]   for i = 0..nf+nc-2  (all temperatures except outer clad, which is prescribed)
//
// Algebraic variables (no time derivative):
//   T[nf+nc-1], qqv, gap, pfc, hgap, all mechanical quantities
// \coderef{axial_layer_state}
struct AxialLayerState {
    int nf, nc;

    // ---- Thermal ----
    std::vector<double> T;     // temperature at each node [K], size nf+nc
    std::vector<double> dT;    // dT/dt [K/s], size nf+nc

    // Power density [W/m3]
    double qqv  = 0.0;

    // ---- Gap ----
    double gap  = 0.0;   // radial gap width [m]
    double pfc  = 0.0;   // pellet-cladding contact pressure [MPa]
    double hgap = 0.0;   // total gap conductance [W/(m2·K)]
    bool   gapOpen = true; // gap state

    // ---- Fuel mechanics ----
    std::vector<double> eft;    // fuel linear thermal expansion strain [-], size nf
    std::vector<double> efh;    // fuel hoop strain [-], size nf
    std::vector<double> efr;    // fuel radial strain [-], size nf
    double              efz;    // fuel axial strain [-] (uniform over layer)
    double              dEfz;   // fuel axial strain rate [1/s]
    std::vector<double> efs;    // fuel volumetric swelling strain [-], size nf (0 for FRED-ROD)
    // Axial strain offset at gas-bond closure: C = efz - ez recorded when gapOpen→closed.
    // Enforced as efz - ez = axialOffset (algebraically equivalent to defz = dez).
    // Unused for liquid bond (axialOffset is always 0 there).
    double              axialOffset = 0.0;

    std::vector<double> sigfh;  // fuel hoop stress [MPa], size nf
    std::vector<double> sigfr;  // fuel radial stress [MPa], size nf
    std::vector<double> sigfz;  // fuel axial stress [MPa], size nf
    std::vector<double> sigf;   // fuel effective stress [MPa], size nf

    // ---- Cladding mechanics ----
    std::vector<double> et;     // clad linear thermal expansion strain [-], size nc
    std::vector<double> eh;     // clad hoop strain [-], size nc
    std::vector<double> er;     // clad radial strain [-], size nc
    double              ez;     // clad axial strain [-]
    double              dEz;    // clad axial strain rate [1/s]
    std::vector<double> ec;     // clad hoop creep strain [-], size nc (0 for FRED-ROD)
    std::vector<double> sigh;   // clad hoop stress [MPa], size nc
    std::vector<double> sigr;   // clad radial stress [MPa], size nc
    std::vector<double> sigz;   // clad axial stress [MPa], size nc
    std::vector<double> sig;    // clad effective stress [MPa], size nc

    // ---- Fuel thermal conductivity (two paths, mutually exclusive) -----------
    //
    // HeatConduction::computeResiduals dispatches on whether k_fuel_per_node
    // is populated:
    //
    //   PATH A — scalar (FRED-ROD, FRED-OX, default):
    //     k_fuel_per_node is empty.
    //     kf = thermalConductivity(T_half) * k_irr_factor
    //     T_half = 0.5*(T[i]+T[i+1]) is the CURRENT Newton iterate, so kf
    //     updates every residual evaluation and the Jacobian sees the live
    //     temperature dependence of conductivity.
    //     k_irr_factor is a scalar irradiation-correction ratio (≤1),
    //     updated once per accepted step in afterAcceptedStep.
    //
    //   PATH B — per-node (FRED-M-Na with Zr redistribution):
    //     k_fuel_per_node[i] holds the absolute irradiated conductivity
    //     [W/(m·K)] at node i, computed from the local post-redistribution
    //     composition (zr_wf, pu_wf) and porosity state at node i.
    //     kf = 0.5*(k_fuel_per_node[i] + k_fuel_per_node[i+1]) for interval i.
    //     IMPORTANT: k_fuel_per_node is filled once in afterAcceptedStep and is
    //     FROZEN for the entire Newton solve of the next step — it does NOT
    //     update with T during Newton iterations (one-step lag on T).  This is
    //     acceptable because the primary benefit is capturing the radial
    //     composition gradient k(zr(r)), not the temperature dependence; the
    //     latter would require re-evaluating the full irradiation model every
    //     residual call, which would need FredMNaNodeState inside HeatConduction.
    //
    // Apps that do not set k_fuel_per_node (FRED-ROD, FRED-OX) use Path A
    // automatically; no code change is required in those apps.
    double              k_irr_factor    = 1.0;  // Path A scalar correction [-]
    std::vector<double> k_fuel_per_node;         // Path B per-node k [W/(m·K)], size nf

    // ---- Deformed geometry ----
    std::vector<double> rad;    // current node radii [m], size nf+nc
    std::vector<double> drf;    // current fuel ring thickness [m], size nf-1 (between nodes)
    std::vector<double> drc;    // current clad ring thickness [m], size nc-1
    double rfi, rfo, rci, rco;  // current key radii [m]

    explicit AxialLayerState(int nf_, int nc_)
        : nf(nf_), nc(nc_),
          T(nf_+nc_, 0.0), dT(nf_+nc_, 0.0),
                    eft(nf_, 0.0), efh(nf_, 0.0), efr(nf_, 0.0), efz(0.0), dEfz(0.0),
                    efs(nf_, 0.0),
          sigfh(nf_, 0.0), sigfr(nf_, 0.0), sigfz(nf_, 0.0), sigf(nf_, 0.0),
                    et(nc_, 0.0), eh(nc_, 0.0), er(nc_, 0.0), ez(0.0), dEz(0.0),
          ec(nc_, 0.0),
          sigh(nc_, 0.0), sigr(nc_, 0.0), sigz(nc_, 0.0), sig(nc_, 0.0),
          rad(nf_+nc_, 0.0), drf(nf_-1, 0.0), drc(nc_-1, 0.0),
          rfi(0.0), rfo(0.0), rci(0.0), rco(0.0)
    {}

    // Update deformed geometry from current strain state.
    // Returns false if any dimension goes negative (solver should reduce time step).
    // \coderef{update_geometry}
    bool updateGeometry(const std::vector<double>& rad0,
                        double drf0, double drc0, double dz0,
                        double& dzf, double& dzc) const
    {
        auto& self = *const_cast<AxialLayerState*>(this);

        // Fuel node radii: rad[i] = rad0[i] * (1 + efh[i])
        // rad[0] = 0 is valid for a solid pellet (centre node); use strict < 0.
        for (int i = 0; i < nf; ++i) {
            self.rad[i] = rad0[i] * (1.0 + efh[i]);
            if (self.rad[i] < 0.0) return false;
        }
        // Fuel ring thicknesses: drf[i] = drf0 * (1 + 0.5*(efr[i]+efr[i+1]))
        for (int i = 0; i < nf-1; ++i) {
            double efr_avg = 0.5 * (efr[i] + efr[i+1]);
            self.drf[i] = drf0 * (1.0 + efr_avg);
            if (self.drf[i] <= 0.0) return false;
        }
        // Fuel axial stretch
        dzf = dz0 * (1.0 + efz);
        if (dzf <= 0.0) return false;

        // Boundary nodes sit exactly on the surfaces (vertex-centred scheme),
        // so surface radii are just the deformed boundary-node radii.
        self.rfi = rad0[0]      * (1.0 + efh[0]);      // 0 for solid pellet
        self.rfo = rad0[nf - 1] * (1.0 + efh[nf - 1]); // outer fuel surface

        // Clad node radii
        for (int i = 0; i < nc; ++i) {
            self.rad[nf+i] = rad0[nf+i] * (1.0 + eh[i]);
            if (self.rad[nf+i] <= 0.0) return false;
        }
        // Clad ring thicknesses
        for (int i = 0; i < nc-1; ++i) {
            double er_avg = 0.5 * (er[i] + er[i+1]);
            self.drc[i] = drc0 * (1.0 + er_avg);
            if (self.drc[i] <= 0.0) return false;
        }
        dzc = dz0 * (1.0 + ez);
        if (dzc <= 0.0) return false;

        self.rci = rad0[nf]          * (1.0 + eh[0]);      // inner clad surface
        self.rco = rad0[nf + nc - 1] * (1.0 + eh[nc - 1]); // outer clad surface
        return true;
    }
};

// Full rod state: one AxialLayerState per axial layer.
// \coderef{rod_state}
struct RodState {
    int nf, nc, nz;
    double gpres;        // internal gas pressure [MPa]
    double mu0;          // initial fill-gas amount [mol]
    double fgrel;        // released fission gas [mol] (0 for FRED-ROD)
    std::vector<AxialLayerState> layers;

    RodState(int nf_, int nc_, int nz_)
        : nf(nf_), nc(nc_), nz(nz_),
          gpres(0.0), mu0(0.0), fgrel(0.0),
          layers(nz_, AxialLayerState(nf_, nc_))
    {}
};

} // namespace fred
