#pragma once
#include "platform/GapMaterial.hpp"
#include <algorithm>
#include <cmath>

namespace fred {

// Sodium bond gap material for FRED-M-Na (U-Pu-Zr metallic fuel pins).
//
// Two modes correspond to the igap flag in Gaphtc.for (Timpano/EPFL FRED_M_OCT24):
//
//   TDependent (igap=1):
//     Temperature-dependent Na thermal conductivity (Touloukian polynomial):
//       k_Na(T) = 93 - 5.81e-2*T_C + 1.173e-5*T_C^2   [W/(m·K)]
//     where T_C = T_K - 273.15.  Gap conductance is k_Na / gap_eff.
//
//   Constant (igap=2):
//     k_Na = 62.9 W/(m·K) (value recommended by ANL for SAS modelling).
//     Conductance is clamped to [1e5, 1e6] W/(m²·K).
//     When gap is soft or closed the conductance is forced to 1e6 W/(m²·K).
//
// isGasBond() returns false: the platform skips the gas-jump model.
// computeMNaGapConductance in FredMNaResiduals uses gapConductivity(T) and
// then calls clampConductance() for the igap=2 bounds.
//
enum class SodiumMode {
    TDependent = 1,  // igap=1: T-dependent polynomial
    Constant   = 2   // igap=2: 62.9 W/(m·K) with conductance caps
};

class Sodium : public GapMaterial {
public:
    explicit Sodium(SodiumMode mode = SodiumMode::TDependent) : m_mode(mode) {}

    // Thermal conductivity of liquid sodium [W/(m·K)]
    // TDependent: Touloukian polynomial (Gaphtc.for, igap=1)
    // Constant  : 62.9 W/(m·K) (ANL SAS recommendation)
    double gapConductivity(double T_K) const override {
        if (m_mode == SodiumMode::Constant) return 62.9;
        const double T_C = std::max(T_K - 273.15, 0.0);
        return 93.0 - 5.81e-2 * T_C + 1.173e-5 * T_C * T_C;
    }

    // Liquid bond — platform must NOT apply the gas-jump model.
    bool isGasBond() const override { return false; }

    // igap=2: clamp h_gap to [1e5, 1e6] W/(m²·K); force 1e6 when closed.
    // igap=1: identity (no clamping).
    double clampConductance(double h_gap, bool gap_closed) const override {
        if (m_mode == SodiumMode::Constant) {
            h_gap = std::max(h_gap, 1.0e5);
            if (h_gap >= 1.0e6 || gap_closed)
                h_gap = 1.0e6;
        }
        return h_gap;
    }

    SodiumMode mode() const { return m_mode; }

private:
    SodiumMode m_mode;
};

} // namespace fred
