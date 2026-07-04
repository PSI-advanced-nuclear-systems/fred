#include "DummyFuelPellet.hpp"
#include "platform/Constants.hpp"

namespace fred {

double DummyFuelPellet::thermalConductivity(double) const  { return 10.0; }
double DummyFuelPellet::heatCapacity(double) const         { return 100.0; }
double DummyFuelPellet::thermalExpansionStrain(double T) const { return 1.0e-5 * (T - T_REF); }
double DummyFuelPellet::youngsModulus(double, double) const { return 100000.0; }
double DummyFuelPellet::poissonRatio() const                { return 0.3; }
double DummyFuelPellet::referenceDensity() const            { return 10000.0; }
double DummyFuelPellet::theoreticalDensity() const          { return 10000.0; }
double DummyFuelPellet::meltingTemperature() const          { return 3000.0; }

} // namespace fred
