#pragma once
#include "FuelPelletMaterial.hpp"
#include "CladdingMaterial.hpp"

namespace fred {

// Gap thermal conductance between fuel outer surface and cladding inner surface.
// Decomposes into three parallel components (Ross & Stoute model):
//
//   h_total = h_gap + h_rad + h_contact
//
// h_gap    : gap-medium conduction term (Lanning-Hann accommodation model)
// h_rad    : linearised radiative transfer between concentric cylinders (only when isGasBond=True)
// h_contact: solid contact conductance (non-zero only when gap is closed)
struct GapConductanceResult {
    double h_gap     = 0.0;  // [W/(m2·K)] medium conduction term
    double h_rad     = 0.0;
    double h_contact = 0.0;
    double h_total   = 0.0;
};

// Gap-medium conductance (Lanning-Hann jump + roughness-limited gap).
// k_gap is the effective thermal conductivity of the medium in the gap
// [W/(m·K)] and is supplied by the application/material model.
// reloc: fuel relocation fraction [-], shrinks the effective gap.
double computeGapMediumConductance(
    double T_fuel, double T_clad,
    double gap_width, double gpres,
    double k_gap,
    double ruff, double rufc,
    double reloc = 0.0);

// Linearised radiation conductance between concentric cylinders.
double computeRadiationConductance(
    double T_fuel, double T_clad,
    double rfo0, double rci0);

// Open-gap conductance helper: medium conduction + optional radiation.
// Kept as a convenience wrapper.
GapConductanceResult computeOpenGapConductance(
    double T_fuel, double T_clad,
    double gap_width, double gpres,
    double k_gap,
    double rfo0, double rci0,
    double ruff, double rufc,
    double reloc = 0.0);

// Full gap conductance (open-gap terms + contact conductance).
GapConductanceResult computeGapConductance(
    double T_fuel, double T_clad,
    double gap_width, double pfc, double gpres,
    double k_gap,
    double rfo0, double rci0,
    double ruff, double rufc,
    const FuelPelletMaterial& fuel, const CladdingMaterial& clad,
    double reloc = 0.0);

// Solid contact conductance [W/(m²·K)] using the Ross & Stoute model.
// Non-zero only when pfc > 0 and gap_width <= sqrt(ruff²+rufc²).
double computeContactConductance(
    double T_fuel, double T_clad,
    double gap_width, double pfc, double gpres,
    double ruff, double rufc,
    const FuelPelletMaterial& fuel, const CladdingMaterial& clad);

} // namespace fred
