#pragma once

namespace fred {

// Abstract interface for gap-medium properties.
//
// All bond types expose a single thermal conductivity k [W/(m·K)].
// The solver converts k to a gap conductance h [W/(m²·K)] using the
// model selected by isGasBond():
//
//   Gas bond  (isGasBond() = true, default):
//     h = computeGapMediumConductance(k, gap, ...) + h_rad
//     Lanning-Hann gas-jump model with actual gap width, plus linearised
//     radiation between concentric cylinders.
//
//   Liquid bond  (isGasBond() = false):
//     h = k / gap_eff   (no gas-jump, no radiation)
//     clampConductance() may be overridden to enforce min/max bounds
//     (e.g. sodium forced to 1e6 W/m²K on gap closure).
//
// \coderef{gap_material_iface}
class GapMaterial {
public:
    virtual ~GapMaterial() = default;

    // Gap-medium thermal conductivity [W/(m·K)] at average gap temperature T [K].
    virtual double gapConductivity(double T) const = 0;

    // Bond type selector.
    // true  = fill-gas (He/Xe/Kr): platform applies Lanning-Hann + radiation.
    // false = liquid bond (Na/Pb/…): platform uses k / gap_eff directly.
    // \coderef{gap_bond_selector}
    virtual bool isGasBond() const { return true; }

    // Optional post-calculation clamp on h_gap [W/(m²·K)].
    // Called by the solver after computing h = k/gap_eff, before storing.
    // Override for liquid-bond materials that impose min/max conductance bounds.
    // Default is identity (no clamping).
    virtual double clampConductance(double h_gap, bool /*gap_closed*/) const {
        return h_gap;
    }

    // See FuelPelletMaterial::isThreadSafe() for the full rationale.
    virtual bool isThreadSafe() const { return true; }
};

} // namespace fred
