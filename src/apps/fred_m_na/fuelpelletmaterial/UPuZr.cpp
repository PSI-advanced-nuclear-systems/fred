#include "UPuZr.hpp"
#include <algorithm>
#include <cmath>

// U-Pu-Zr metallic fuel property correlations.
// Sources: Flamb.for / Ftexp.for / Felmod.for / Fpoir.for / Fcp.for /
//          Fcreep.for (upuzr branch) from Timpano/EPFL FRED_M_OCT24.SRC.
//          Karahan (2009), Hofman (1985), Gruber.

namespace {

using std::max; using std::min; using std::exp; using std::pow; using std::sqrt;

// Thermal conductivity of fresh U-Pu-Zr [W/(m·K)]  — Flamb.for, Aydin/Karahan
static double upuzrFreshConductivity(double T_K, double pu, double zr) {
    const double aak = 17.5 * ((1.0 - 2.23 * zr) / (1.0 + 1.61 * zr) - 2.62 * pu);
    const double bbk = 1.54e-2 * (1.0 + 0.061 * zr) / (1.0 + 1.61 * zr);
    const double cck = 9.38e-6 * (1.0 - 2.7 * pu);
    return aak + bbk * T_K + cck * T_K * T_K;
}

// Irradiation-corrected thermal conductivity [W/(m·K)]  — Flamb.for, f=1
// Accounts for sodium infiltration (Maxwell-Eucken) and gas porosity.
// psod : fraction due to sodium infiltration. 
static double thermalConductivitySodiumInfiltration(double T_K, double pu, double zr,
                                                double poros_tot, double poros_gas,
                                                double psod)
{
    const double ulmb = upuzrFreshConductivity(T_K, pu, zr);

    // Sodium thermal conductivity [W/(m·K)]  — MFUEL correlation
    const double nalamb = 110.45 - 0.065112 * T_K + 0.00001543 * T_K * T_K
                        - 0.0000000024617 * T_K * T_K * T_K;

    const double pna  = psod * (poros_tot - poros_gas);
    const double pgas = poros_gas + (poros_tot - poros_gas) * (1.0 - psod);

    const double na_u = nalamb / max(ulmb, 1.0e-6);
    const double pc   = 1.0 - 3.0 * pna / max(1.0 - pgas, 1.0e-6)
                      * (1.0 - na_u) / (1.163 + 1.837 * na_u);
    return ulmb * pc * pow(max(1.0 - pgas, 0.0), 1.5);
}

// Empirical burnup-correction thermal conductivity [W/(m·K)]  — Flamb.for, f=2
static double thermalConductivityEmpirical(double T_K, double pu, double zr,
                                                  double bup_FIMA)
{
    const double ulmb = upuzrFreshConductivity(T_K, pu, zr);
    const double pp = 0.135 * bup_FIMA * 100.0;
    double correction;
    if (bup_FIMA <= 0.02)
        correction = (1.0 - pp) / (1.0 + 1.7 * pp);
    else if (bup_FIMA <= 0.05)
        correction = 0.5 + 0.0667 * (bup_FIMA * 100.0 - 2.0);
    else
        correction = 0.7;
    return ulmb * correction;
}

// ESFR-SIMPLE sigmoid-fit correction  — Flamb.for, f=3
// Fitted correction factor as difference of two logistic curves in burnup (at%).
static double thermalConductivityEsfrSimple(double T_K, double pu, double zr,
                                                   double bup_FIMA)
{
    const double ulmb = upuzrFreshConductivity(T_K, pu, zr);
    const double bup_atpct = bup_FIMA * 100.0;
    const double a1 = 100.27110592, b1 = 0.5474772,  c1 = -8.15852445;
    const double a2 =  99.62701,    b2 = 0.64761187, c2 = -6.4644582;
    const double correction = a1 / (1.0 + exp(-b1 * (bup_atpct - c1)))
                            - a2 / (1.0 + exp(-b2 * (bup_atpct - c2)));
    return ulmb * correction;
}

// Linear thermal expansion strain relative to 298.15 K [-]  — Ftexp.for (Karahan 2009)
static double upuzrThermalExpansionStrain(double T_K) {
    if (T_K <= 877.0)
        return 1.76e-5 * (T_K - 298.15);
    else if (T_K >= 936.0)
        return 2.01e-5 * (T_K - 298.15);
    else {
        const double slope = (2.01e-5 - 1.76e-5) / (936.0 - 877.0);
        return (slope * (T_K - 877.0) + 1.76e-5) * (T_K - 298.15);
    }
}

// Young's modulus [Pa]  — Felmod.for (Karahan 2009)
static double upuzrElasticModulus(double T_K, double density, double rho0,
                                   double pu, double zr)
{
    const double por = 1.0 - density / rho0;
    double E = 1.6e11;
    E *= (1.0 - 1.2 * por);
    E *= ((1.0 + 0.17 * zr) / (1.0 + 1.34 * zr) - pu);
    E *= (1.0 - 1.06 * (T_K - 588.0) / 1405.0);
    return max(E, 1.0e8);
}

// Poisson's ratio [-]  — Fpoir.for (Karahan 2009)
static double upuzrPoissonRatio(double T_K, double density, double rho0,
                                 double /*pu*/, double zr)
{
    const double por = 1.0 - density / rho0;
    double nu = 0.24;
    nu *= (1.0 - 0.8 * por);
    nu *= (1.0 + 3.4 * zr) / (1.0 + 1.9 * zr);
    nu *= (1.0 + 1.2 * (T_K - 588.0) / 1405.0);
    return max(0.01, min(0.49, nu));
}

// Heat capacity [J/(kg·K)]  — Fcp.for (Karahan 2009)
static double upuzrHeatCapacity(double T_K, double pu, double zr) {
    const double tc = T_K - 273.15;
    const double ma = (238.02891 * (1.0 - pu - zr) + 244.06 * pu + 91.22 * zr) * 1.0e-3;
    double cp;
    if (tc <= 600.0)
        cp = 26.58 / ma + 0.027 * tc / ma;
    else if (tc >= 650.0)
        cp = 15.84 / ma + 0.026 * tc / ma;
    else {
        const double cp1 = 26.58 / ma + 0.027 * 600.0 / ma;
        const double cp2 = 15.84 / ma + 0.026 * 650.0 / ma;
        cp = (cp2 - cp1) / 50.0 * (tc - 600.0) + cp1;
    }
    return cp;
}

// Effective creep strain rate [1/s]  — Fcreep.for, upuzr branch (Karahan 2009)
//   sig_eff    : von Mises effective stress [Pa]
//   sgh,sgz,sgr: stress components [Pa]
//   gasp_Pa    : gas pressure [Pa]
//   dotfissden : fission rate density [fiss/m3/s]
//   swopen     : open-pore swelling (GRSIS; 0 → simplified)
static double upuzrCreepRate(double T_K, double sig_eff,
                              double sgh, double sgz, double sgr,
                              double gasp_Pa, double dotfissden, double swopen)
{
    if (sig_eff <= 0.0) return 0.0;

    const double shyd = max(sgh + sgz + sgr + 3.0 * gasp_Pa, 0.0);
    const double por  = max(swopen, 0.0);
    double ac = 0.0;
    if (por > 0.0)
        ac = (1.0 / 6.0) * pow(min(por, 0.1) / 0.1, 1.5);

    const double seq = sqrt(0.5 * ((sgh - sgz)*(sgh - sgz)
                                  + (sgh - sgr)*(sgh - sgr)
                                  + (sgz - sgr)*(sgz - sgr))
                            + 3.0 * ac * shyd * shyd);
    const double s_MPa = seq * 1.0e-6;

    double fcreep;
    if (T_K < 923.15) {
        fcreep = (5.0e3 * s_MPa + 6.0 * pow(s_MPa, 4.5)) * exp(-26170.0 / T_K)
               + 7.7e-23 * dotfissden * 1.0e-6 * s_MPa;
    } else {
        fcreep = 0.08 * pow(s_MPa, 3.0) * exp(-14350.0 / T_K);
    }
    return max(fcreep, 0.0);
}

// Solidus temperature [K] for U-Pu-Zr — Mmelt.for (Karahan 2009)
// 25-point table at Pu=13.25 wt%, Zr from 5–29 wt%.
// Nearest-neighbour lookup on Zr weight percent.
static double upuzrSolidusTemperature(double zr_wf) {
    static const double zrwtab[25] = {
         5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29
    };
    static const double tsoltab[25] = {
        1382.22,1394.33,1404.02,1413.70,1423.39,1431.87,1441.55,1451.24,
        1462.14,1471.83,1483.94,1496.05,1508.16,1520.27,1534.80,1546.91,
        1561.44,1575.97,1589.29,1602.61,1617.14,1630.46,1643.78,1657.10,
        1670.42
    };
    const double zr_wtpct = zr_wf * 100.0;
    int best = 0;
    double best_diff = std::abs(zr_wtpct - zrwtab[0]);
    for (int i = 1; i < 25; ++i) {
        const double d = std::abs(zr_wtpct - zrwtab[i]);
        if (d < best_diff) { best_diff = d; best = i; }
    }
    return tsoltab[best];
}

} // anonymous namespace

// ---------------------------------------------------------------------------
namespace fred {

UPuZr::UPuZr(double pu_weight_frac, double zr_weight_frac, double rho0)
    : m_pu(pu_weight_frac), m_zr(zr_weight_frac)
{
    const double ur = 1.0 - m_pu - m_zr;
    m_tden = 1.0 / (ur / 19040.0 + m_pu / 19800.0 + m_zr / 6510.0);
    m_rho0 = (rho0 > 0.0) ? rho0 : 0.75 * m_tden;
}

double UPuZr::thermalConductivity(double T) const {
    return upuzrFreshConductivity(T, m_pu, m_zr);
}

double UPuZr::thermalConductivityIrradiated(double T_K, double bup_FIMA,
                                              double poros_tot, double poros_gas,
                                              double psod) const
{
    switch (m_conductivity_model) {
        case ConductivityModel::EmpiricalBurnup:
            return thermalConductivityEmpirical(T_K, m_pu, m_zr, bup_FIMA);
        case ConductivityModel::EsfrSimple:
            return thermalConductivityEsfrSimple(T_K, m_pu, m_zr, bup_FIMA);
        default: // DetailedNaSodium
            return thermalConductivitySodiumInfiltration(T_K, m_pu, m_zr, poros_tot, poros_gas, psod);
    }
}

double UPuZr::heatCapacity(double T) const {
    return upuzrHeatCapacity(T, m_pu, m_zr);
}

double UPuZr::thermalExpansionStrain(double T) const {
    return upuzrThermalExpansionStrain(T);
}

double UPuZr::youngsModulus(double T, double density) const {
    return upuzrElasticModulus(T, density, m_rho0, m_pu, m_zr) * 1.0e-6; // Pa → MPa
}

// Returns baseline value at reference conditions; use poissonRatioFull() for T-dependence.
double UPuZr::poissonRatio() const {
    return 0.24;
}

double UPuZr::referenceDensity()   const { return m_rho0; }
double UPuZr::theoreticalDensity() const { return m_tden; }

double UPuZr::meltingTemperature() const {
    return upuzrSolidusTemperature(m_zr);
}

double UPuZr::solidusTemperature(double zr_wf) const {
    return upuzrSolidusTemperature(zr_wf);
}

// UPuZr-specific extended methods (beyond FuelPelletMaterial interface).

double UPuZr::poissonRatioFull(double T_K, double density) const {
    return upuzrPoissonRatio(T_K, density, m_rho0, m_pu, m_zr);
}

double UPuZr::creepRate(double T_K, double sig_eff,
                         double sgh, double sgz, double sgr,
                         double gasp_Pa, double dotfissden, double swopen) const
{
    return upuzrCreepRate(T_K, sig_eff, sgh, sgz, sgr, gasp_Pa, dotfissden, swopen);
}

} // namespace fred
