#pragma once

namespace fred {

// Abstract interface for fuel pellet material properties.
//
// Only properties that are meaningful without irradiation are defined here.
// Density- and burnup-dependent variants (swelling, fission-gas conductivity
// correction, etc.) are introduced by FredOxFuelPelletMaterial in FRED-OX,
// which inherits from this class and adds the irradiation-specific overrides.
//
// See doc/architecture.md for the full inheritance hierarchy.
// \coderef{fuel_pellet_iface}
class FuelPelletMaterial {
public:
    virtual ~FuelPelletMaterial() = default;

    // Thermal conductivity [W/(m·K)] at temperature T [K].
    // As-fabricated (unirradiated) correlation; density and burnup dependence
    // is introduced in FredOxFuelPelletMaterial.
    virtual double thermalConductivity(double T) const = 0;

    // Specific heat capacity [J/(kg·K)]
    virtual double heatCapacity(double T) const = 0;

    // Linear thermal expansion strain relative to 293.15 K [-]
    virtual double thermalExpansionStrain(double T) const = 0;

    // Young's modulus [MPa] — density dependence retained (porosity correction)
    virtual double youngsModulus(double T, double density) const = 0;

    // Poisson's ratio [-]
    virtual double poissonRatio() const = 0;

    // As-fabricated density [kg/m3]
    virtual double referenceDensity() const = 0;

    // Theoretical (fully dense) density [kg/m3]
    virtual double theoreticalDensity() const = 0;

    // Melting temperature [K]
    virtual double meltingTemperature() const = 0;

    // True if this instance's virtual methods may be safely called
    // concurrently from multiple OpenMP worker threads (one per axial
    // layer). Default true: all built-in C++ material implementations are
    // pure functions of their arguments with no shared mutable state.
    // pybind11 trampolines for Python-subclassed materials override this to
    // false, since calling back into the Python interpreter from a
    // non-GIL-holding worker thread would corrupt or crash it. RodResiduals
    // clamps threads back to 1 whenever any bound material reports false.
    virtual bool isThreadSafe() const { return true; }
};

} // namespace fred
