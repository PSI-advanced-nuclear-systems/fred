#pragma once
#include "platform/GapMaterial.hpp"
#include "platform/Constants.hpp"
#include <algorithm>
#include <cmath>

namespace fred {

// Dummy gap material for verification tests.
//
// Uses He conductivity for the Lanning-Hann gap conductance model
// (same as He() but bundled as a single self-contained class).
// isGasBond() = true: (1) gap annular volume counted in as-fabricated free-gas volume, (2) radiative conductance term h_rad is enabled.
class DummyGapMaterial : public GapMaterial {
public:
    double gapConductivity(double T) const override {
        return 2.639e-3 * std::pow(T, 0.7085);
    }

    bool isGasBond() const override { return true; }
};

} // namespace fred
