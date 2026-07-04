#pragma once

namespace fred {

// Abstract interface for cladding material properties.
//
// Covers thermo-elastic properties without irradiation effects.
// Irradiation-induced creep, void swelling, and embrittlement are
// introduced by derived classes in FRED-OX and FRED-M-Na.
//
// See doc/architecture.md for the full inheritance hierarchy.
// \coderef{cladding_iface}
class CladdingMaterial {
public:
    virtual ~CladdingMaterial() = default;

    // Thermal conductivity [W/(m·K)]
    virtual double thermalConductivity(double T) const = 0;

    // Specific heat capacity [J/(kg·K)]
    virtual double heatCapacity(double T) const = 0;

    // Linear thermal expansion strain relative to 293.15 K [-]
    virtual double thermalExpansionStrain(double T) const = 0;

    // Young's modulus [MPa]
    virtual double youngsModulus(double T) const = 0;

    // Poisson's ratio [-]
    virtual double poissonRatio() const = 0;

    // Meyer hardness [Pa] — used in solid contact conductance (Ross & Stoute)
    virtual double meyerHardness(double T) const = 0;

    // As-fabricated density [kg/m3]
    virtual double referenceDensity() const = 0;

    // Equivalent hoop creep strain rate [1/s] at temperature T [K] and hoop stress sigma [MPa].
    // Default: zero (FRED-ROD cladding materials — no irradiation creep model).
    // Override in FRED-OX cladding classes.
    // \coderef{cladding_creep_hook}
    virtual double creepRate(double T, double sigma) const { (void)T; (void)sigma; return 0.0; }

    // See FuelPelletMaterial::isThreadSafe() for the full rationale.
    virtual bool isThreadSafe() const { return true; }
};

} // namespace fred
