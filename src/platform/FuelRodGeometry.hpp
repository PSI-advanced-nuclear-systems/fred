#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include "Constants.hpp"

namespace fred {

// Initial (as-fabricated) geometry of a cylindrical fuel rod.
// All lengths in metres.  Indexing throughout is 0-based.
//
// Radial discretisation (vertex-centred, matching the original Fortran FRED):
//   Fuel nodes  : i = 0 .. nf-1   (node 0 at rfi0, node nf-1 at rfo0)
//   Clad nodes  : i = nf .. nf+nc-1 (node nf at rci0, node nf+nc-1 at rco0)
//
// Nodes sit exactly on the physical surfaces, so pressure boundary conditions
// are applied at the correct radii.  nf ≥ 2 and nc ≥ 2 are required.
//
//   drf0 = (rfo0 - rfi0) / (nf - 1)   — spacing between adjacent fuel nodes
//   drc0 = (rco0 - rci0) / (nc - 1)   — spacing between adjacent clad nodes
//   rad0[i]       = rfi0 + i * drf0   (fuel)
//   rad0[nf+i]    = rci0 + i * drc0   (clad)
//   rad_half[i]   = midpoint(rad0[i], rad0[i+1])  (face between nodes i and i+1)
//
// Control-volume areas (half-cell at boundary nodes, full-cell at interior):
//   inner/outer boundary nodes extend only to drf0/2 (or drc0/2) inward/outward.
// \coderef{fuel_rod_geometry}
struct FuelRodGeometry {
    int nf;   // number of fuel radial nodes  (≥ 2)
    int nc;   // number of cladding radial nodes  (≥ 2)
    int nz;   // number of axial layers

    // Key radii [m]
    double rfi0;  // inner fuel radius  (0 for solid pellet)
    double rfo0;  // outer fuel radius
    double rci0;  // inner cladding radius
    double rco0;  // outer cladding radius

    // Derived uniform spacings [m]
    double drf0;  // fuel node spacing    = (rfo0 - rfi0) / (nf - 1)
    double drc0;  // clad node spacing    = (rco0 - rci0) / (nc - 1)

    // Axial layer heights [m], size nz (may be non-uniform)
    std::vector<double> dz0;

    // Node radii [m], size nf+nc
    //   rad0[i]    = rfi0 + i * drf0               for i = 0..nf-1
    //   rad0[nf+i] = rci0 + i * drc0               for i = 0..nc-1
    std::vector<double> rad0;

    // Face (half-node) radii [m], size nf+nc-1
    //   rad_half[i] = (rad0[i] + rad0[i+1]) / 2    for i = 0..nf-2  (fuel-fuel)
    //   rad_half[nf-1]                              skip (gap, handled separately)
    //   rad_half[nf+i] = (rad0[nf+i] + rad0[nf+i+1]) / 2  for i = 0..nc-2  (clad-clad)
    std::vector<double> rad_half;

    // Control-volume cross-sectional areas [m2], size nf+nc
    //   Boundary nodes: half-cell width → pi*((r+dref/2)^2 - r^2) or pi*(r^2-(r-dref/2)^2)
    //   Interior nodes: full-cell width → 2*pi*r*dref
    std::vector<double> area0;

    double azf0;  // total fuel cross-sectional area [m2] = pi*(rfo0^2 - rfi0^2)
    double azc0;  // total cladding cross-sectional area [m2]

    // Surface roughnesses [m] (used in gap conductance)
    double ruff = 1.0e-6;   // fuel outer surface roughness
    double rufc = 1.0e-6;   // cladding inner surface roughness

    // Gas plenum volume [m3]
    double vgp = 0.0;

    // Build all derived arrays from the key inputs.
    // \coderef{geometry_build}
    void build() {
        if (nf < 2 || nc < 2)
            throw std::runtime_error("FuelRodGeometry: nf and nc must each be >= 2");
        const int nr = nf + nc;

        drf0 = (rfo0 - rfi0) / (nf - 1);
        drc0 = (rco0 - rci0) / (nc - 1);

        // Node positions: boundary nodes sit exactly on the surfaces
        rad0.resize(nr);
        for (int i = 0; i < nf; ++i) rad0[i]      = rfi0 + i * drf0;
        for (int i = 0; i < nc; ++i) rad0[nf + i]  = rci0 + i * drc0;

        // Face radii (mid-points between adjacent nodes)
        rad_half.resize(nr - 1);
        for (int i = 0; i < nf - 1; ++i)
            rad_half[i] = 0.5 * (rad0[i] + rad0[i + 1]);
        // rad_half[nf-1] is skipped (gap)
        for (int i = 0; i < nc - 1; ++i)
            rad_half[nf + i] = 0.5 * (rad0[nf + i] + rad0[nf + i + 1]);

        // Control-volume areas
        area0.resize(nr);
        // Fuel — inner boundary half-cell, interior full-cells, outer boundary half-cell
        auto halfCell = [](double rInner, double rOuter) {
            return PI * (rOuter * rOuter - rInner * rInner);
        };
        area0[0]      = halfCell(rfi0, rfi0 + 0.5 * drf0);
        for (int i = 1; i < nf - 1; ++i)
            area0[i] = halfCell(rad0[i] - 0.5 * drf0, rad0[i] + 0.5 * drf0);
        area0[nf - 1] = halfCell(rfo0 - 0.5 * drf0, rfo0);
        // Cladding
        area0[nf]         = halfCell(rci0, rci0 + 0.5 * drc0);
        for (int i = 1; i < nc - 1; ++i)
            area0[nf + i] = halfCell(rad0[nf + i] - 0.5 * drc0, rad0[nf + i] + 0.5 * drc0);
        area0[nf + nc - 1] = halfCell(rco0 - 0.5 * drc0, rco0);

        azf0 = PI * (rfo0 * rfo0 - rfi0 * rfi0);
        azc0 = PI * (rco0 * rco0 - rci0 * rci0);

        if (dz0.empty())
            dz0.assign(nz, 0.0);
    }
};

} // namespace fred
