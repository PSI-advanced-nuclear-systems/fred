#include "GapModel.hpp"
#include "Constants.hpp"
#include <cmath>
#include <algorithm>

namespace fred {

double computeContactConductance(
    double T_fuel, double T_clad,
    double gap_width, double pfc, double gpres,
    double ruff, double rufc,
    const FuelPelletMaterial& fuel, const CladdingMaterial& clad)
{
    double rough = std::sqrt(ruff * ruff + rufc * rufc);
    if (pfc <= 0.0 || gap_width > rough) return 0.0;

    double kf = fuel.thermalConductivity(T_fuel);
    double kc = clad.thermalConductivity(T_clad);
    double k_contact = 2.0 * kf * kc / (kf + kc);
    double H = clad.meyerHardness(T_clad);
    double delta_p = std::max(pfc - gpres, 0.0) * 1.0e6; // net contact pressure [Pa]
    return 33.3 * k_contact * delta_p / (H * std::sqrt(rough));
}

double computeGapMediumConductance(
    double T_fuel, double T_clad,
    double gap_width, double gpres,
    double k_gap,
    double ruff, double rufc,
    double reloc)
{
    double T_avg  = 0.5 * (T_fuel + T_clad);
    double alpha  = std::max(0.01, 0.425 - 2.3e-4 * std::min(T_avg, 1000.0));
    double p_Pa   = std::max(gpres, 0.1) * 1.0e6;
    double g_jump = 0.024688 * k_gap * std::sqrt(T_avg) * 2.0 / (p_Pa * alpha);

    double rough   = std::sqrt(ruff * ruff + rufc * rufc);
    double gap_eff = std::max(gap_width, rough);

    return k_gap / (gap_eff * (1.0 - reloc) + g_jump);
}

double computeRadiationConductance(
    double T_fuel, double T_clad,
    double rfo0, double rci0)
{
    double emissf      = 0.8, emissc = 0.9;
    double view_factor = 1.0 / (1.0/emissf + (rfo0/rci0) * (1.0/emissc - 1.0));
    return SIGMA_SB * view_factor
           * (T_fuel*T_fuel + T_clad*T_clad) * (T_fuel + T_clad);
}

GapConductanceResult computeOpenGapConductance(
    double T_fuel, double T_clad,
    double gap_width, double gpres,
    double k_gap,
    double rfo0, double rci0,
    double ruff, double rufc,
    double reloc)
{
    GapConductanceResult res;
    res.h_gap = computeGapMediumConductance(
        T_fuel, T_clad,
        gap_width, gpres,
        k_gap,
        ruff, rufc,
        reloc);
    res.h_rad = computeRadiationConductance(T_fuel, T_clad, rfo0, rci0);
    res.h_total = res.h_gap + res.h_rad;
    return res;
}

GapConductanceResult computeGapConductance(
    double T_fuel, double T_clad,
    double gap_width, double pfc, double gpres,
    double k_gap,
    double rfo0, double rci0,
    double ruff, double rufc,
    const FuelPelletMaterial& fuel, const CladdingMaterial& clad,
    double reloc)
{
    GapConductanceResult res = computeOpenGapConductance(
        T_fuel, T_clad, gap_width, gpres, k_gap, rfo0, rci0, ruff, rufc,
        reloc);

    res.h_contact = computeContactConductance(
        T_fuel, T_clad, gap_width, pfc, gpres, ruff, rufc, fuel, clad);

    res.h_total = res.h_gap + res.h_rad + res.h_contact;
    return res;
}

} // namespace fred
