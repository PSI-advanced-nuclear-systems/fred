#pragma once
#include "platform/GapMaterial.hpp"
#include "platform/Constants.hpp"
#include <algorithm>
#include <cmath>

namespace fred {

// As-fabricated helium gap material.
//
// Gas thermal conductivity (power-law fit):
//   k_He(T) = 2.639e-3 * T^0.7085   [W/(m·K)]
//   Source: FRED.f90 gaphtc; widely attributed to tabular fits of NIST data.
//
// Gas bond (isGasBond() = true): the solver applies the Lanning-Hann gas-jump
// model to convert k_He to a gap conductance using the actual gap width.
//
class He : public GapMaterial {
public:
    // He thermal conductivity [W/(m·K)] used by shared gap model blocks.
    double gapConductivity(double T) const override {
        return 2.639e-3 * std::pow(T, 0.7085);
    }
};

} // namespace fred
