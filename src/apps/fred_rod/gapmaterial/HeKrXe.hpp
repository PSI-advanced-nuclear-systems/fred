#pragma once
#include "platform/GapMaterial.hpp"
#include <cmath>
#include <algorithm>

namespace fred {

// He-Kr-Xe mixed-gas gap material for base-irradiated gap conditions.
//
// Constructor parameter: local fuel burnup in atomic percent (at%) of heavy
// metal atoms fissioned.  The mole fractions are derived from a molar
// balance — the same physics as FRED-OX — using typical pin parameters:
//
//   Step 1 — He fill-gas inventory (ideal gas law):
//       n_He0 = P_fill * V_free / (R * T_ref)
//       Typical: P_fill=0.1 MPa, V_free=1.5 cm³, T_ref=293.15 K
//               → n_He0 ≈ 6.15e-5 mol
//
//   Step 2 — Fission gas generated and released:
//       fggen = 0.25 * (bup_atpct/100) * m_fuel / M_HM
//               (0.25 gas atoms per fission; m_fuel=10 g, M_HM=238 g/mol)
//       fgr   = Waltar-Reynolds FGR fraction (restructured zone, upper bound):
//               a = 4.7/bup_MWd * (1 - exp(-bup_MWd/5.9)),  fgr = max(0, 1-a)
//       fgrel = fgr * fggen
//
//   Step 3 — Species moles (FRED.f90 gaphtc FG composition):
//       n_He = n_He0 + 0.0385 * fgrel
//       n_Kr =         0.0769 * fgrel
//       n_Xe =         0.8846 * fgrel
//
//   Step 4 — Mole fractions and geometric-mean conductivity:
//       y_i = n_i / (n_He0 + fgrel),   k_mix = k_He^y_He * k_Kr^y_Kr * k_Xe^y_Xe
//
// Characteristic states (typical pin):
//   bup = 0 at%  → pure He (as-fabricated)
//   bup = 1 at%  (≈ 9.4 MWd/kgHM)  → y_Xe ≈ 0.45  (moderate base irrad.)
//   bup = 5 at%  (≈ 47  MWd/kgHM)  → y_Xe ≈ 0.78  (high burnup)
//   bup → ∞      → 88.46 % Xe / 7.69 % Kr / 3.85 % He  (FG-dominated gap)
class HeKrXe : public GapMaterial {
public:
    // bup_atpct: local burnup in atomic percent (at%).
    // Default (1.0 at%) represents a moderate base-irradiation state.
    explicit HeKrXe(double bup_atpct = 1.0) {
        // 1 at% ≈ 9.4 MWd/kgHM (oxide fuel approximation)
        const double bup = bup_atpct * 9.4;

        // Waltar-Reynolds FGR fraction — restructured zone (upper bound).
        double fgr = 0.0;
        if (bup > 0.0) {
            const double a = 4.7 / bup * (1.0 - std::exp(-bup / 5.9));
            fgr = std::max(0.0, 1.0 - a);
        }

        // He fill gas: ideal gas law at fabrication (P=0.1 MPa, V=1.5 cm³, T=293.15 K)
        const double n_He0 = 0.1e6 * 1.5e-6 / (8.314 * 293.15);  // ~6.15e-5 mol

        // Fission gas generated: 0.25 atoms/fission, 10 g fuel, M_HM = 238 g/mol
        const double fggen   = 0.25 * (bup_atpct / 100.0) * 0.010 / 0.238;
        const double fgrel   = fgr * fggen;

        // FG composition: 88.46% Xe / 7.69% Kr / 3.85% He (FRED.f90 gaphtc)
        const double n_Xe    = 0.8846 * fgrel;
        const double n_Kr    = 0.0769 * fgrel;
        const double n_He_FG = 0.0385 * fgrel;

        const double n_total = n_He0 + n_He_FG + n_Kr + n_Xe;
        m_y_He = (n_He0 + n_He_FG) / n_total;
        m_y_Kr = n_Kr / n_total;
        m_y_Xe = n_Xe / n_total;
    }

    double gapConductivity(double T) const override {
        const double k_He = 2.639e-3 * std::pow(T, 0.7085);
        const double k_Kr = 8.247e-5 * std::pow(T, 0.8363);
        const double k_Xe = 4.351e-5 * std::pow(T, 0.8618);
        return std::pow(k_He, m_y_He) * std::pow(k_Kr, m_y_Kr) * std::pow(k_Xe, m_y_Xe);
    }

private:
    double m_y_He;
    double m_y_Kr;
    double m_y_Xe;
};

} // namespace fred
