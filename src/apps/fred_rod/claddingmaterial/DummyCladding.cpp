#include "DummyCladding.hpp"
#include "platform/Constants.hpp"

namespace fred {

double DummyCladding::thermalConductivity(double) const      { return 10.0; }
double DummyCladding::heatCapacity(double) const             { return 100.0; }
double DummyCladding::thermalExpansionStrain(double T) const { return 1.0e-5 * (T - T_REF); }
double DummyCladding::youngsModulus(double) const            { return 100000.0; }
double DummyCladding::poissonRatio() const                   { return 0.3; }
double DummyCladding::meyerHardness(double) const            { return 1.0e9; }
double DummyCladding::referenceDensity() const               { return 10000.0; }

} // namespace fred
