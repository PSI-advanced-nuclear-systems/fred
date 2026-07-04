#pragma once
#include "platform/SubchannelCoolantProperties.hpp"
#include <cmath>

namespace fred {

// Sodium (Na) thermophysical properties for subchannel cooling.
// Correlations from Baseir.for (FRED-M Timpano reference), lines ~270-287.
class FredMNaSodiumProperties : public SubchannelCoolantProperties {
public:
    // Specific heat [J/(kg*K)] — constant for liquid Na (Baseir.for, Na block)
    double cp(double /*T*/) const override { return 1270.0; }

    // Density [kg/m3] — Baseir.for Na block
    double rho(double T) const override {
        const double f = 1.0 - T / 2503.7;
        return 219.0 + 275.32 * f + 511.58 * std::sqrt(f);
    }

    // Thermal conductivity [W/(m*K)] — Baseir.for Na block
    double k(double T) const override {
        return 1.1641922e2 * std::exp(-7.0184779e-4 * T);
    }

    // Dynamic viscosity [Pa*s] — not used in Peclet-number HTC correlations for Na
    // TODO: add correlation if needed for non-Peclet HTC models
    double mu(double /*T*/) const override { return 0.0; }
};

} // namespace fred
