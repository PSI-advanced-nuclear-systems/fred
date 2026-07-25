#include "HT9.hpp"
#include <cmath>
#include <algorithm>

// HT-9 ferritic-martensitic steel property correlations.
// Sources: Clamb.for / Ctexp.for / Celmod.for / Cpoir.for / Ccreep.for /
//          Cswel.for / Csigy.for from Timpano/EPFL FRED_M_OCT24.SRC.
//          Karahan (2007), Hofmann (1985).
//          D.1.4 ESFR-SIMPLE metallic fuel study: performance and safety analysis


namespace {

using std::min; using std::max; using std::exp; using std::pow; using std::log;

// Thermal conductivity [W/(m·K)]  — Clamb.for
// Hofman 1985 thermal conductivity 
double ThermalConductivity(double T_K) {
    return 29.65 - 6.668e-2 * T_K + 2.184e-4 * T_K * T_K
         - 2.527e-7 * T_K * T_K * T_K + 9.621e-11 * T_K * T_K * T_K * T_K;
}

// Linear thermal expansion strain relative to 293.15 K [-]  — Ctexp.for (Karahan 2007)
double ThermalExpansionStrain(double T_K) {
    return 1.0e-2 * (-0.2191 + 5.678e-4 * T_K + 8.111e-7 * T_K * T_K
                     - 2.576e-10 * T_K * T_K * T_K);
}

// Young's modulus [Pa]  — Celmod.for
double ElasticModulus(double T_K) {
    const double tc = min(T_K - 273.15, 800.0);
    return 2.137e11 - 1.0274e8 * (tc + 273.15);
}

// Poisson's ratio [-]  — Cpoir.for
// Note: the base CladdingMaterial::poissonRatio() takes no T argument.
// This helper is provided for direct use where temperature-dependence matters.
double PoissonRatio(double T_K) {
    const double tc = min(T_K - 273.15, 800.0);
    return 0.5 * (2.137e5 - 102.74 * tc) / (8.964e4 - 53.78 * tc) - 1.0;
}

// Specific heat capacity [J/(kg·K)] 
// Yamanouchi 1992 
double HeatCapacity(double T_K) {
    if (T_K < 800.15) return 416.642 + 0.167*T_K;
    return 69.910 + 0.600*T_K; 
}

// Yield stress [Pa]  — Csigy.for
double YieldStress(double T_K) {
    static const double temht9[] = {293.15, 673.15, 773.15, 873.15, 973.15,
                                     1073.15, 1173.15, 1273.15, 1683.15, 1708.15};
    static const double syht9[]  = {6.69e8, 6.40e8, 5.22e8, 4.885e8, 3.69e8,
                                     2.48e8, 7.95e7, 4.80e7, 1.0e6, 1.0e5};
    constexpr int n = 10;
    if (T_K <= temht9[0])   return syht9[0];
    if (T_K >= temht9[n-1]) return syht9[n-1];
    for (int i = 0; i < n - 1; ++i) {
        if (T_K <= temht9[i + 1]) {
            const double f = (T_K - temht9[i]) / (temht9[i + 1] - temht9[i]);
            return syht9[i] + f * (syht9[i + 1] - syht9[i]);
        }
    }
    return syht9[n - 1];
}

// Meyer hardness [Pa]  — Tabor relation: H ~ 3 * sigma_y
double MeyerHardness(double T_K) {
    return 3.0 * YieldStress(T_K);
}

// Creep strain rate [1/s]  — Ccreep.for
//   sig_Pa  : effective stress [Pa]
//   qqv     : volumetric power density [W/m3] (neutron flux proxy)
//   time_s  : elapsed time [s]
double CreepRate(double T_K, double sig_Pa, double qqv, double time_s) {
    const double sg = sig_Pa * 1.0e-6; // MPa
    // Neutron flux proxy: nflux = qqv * 8.4e5 * 1e-22  [10^22 n/cm2/s]
    const double nflux = qqv * 8.4e5 * 1.0e-22;
    // Irradiation creep (Ccreep.for)
    const double ircreep = (1.83e-4 + 2.59e14 * exp(-7.3e4 / (1.987 * T_K)))
                          * nflux * pow(sg, 1.3);
    // Thermal creep: primary + secondary + tertiary
    const double tem1 = 13.4     * exp(-15027.0  / (1.987 * T_K)) * sg;
    const double tem2 = 8.34e-3  * exp(-26451.0  / (1.987 * T_K)) * pow(sg, 4.0);
    const double tem3 = 4.08e18  * exp(-89167.0  / (1.987 * T_K)) * pow(sg, 0.5);
    const double tem4 = 1.6e-6   * exp(-1.6e-6 * time_s);
    const double tem5 = 1.17e9   * exp(-83142.0  / (1.987 * T_K)) * sg * sg;
    const double tem6 = 8.33e9   * exp(-108276.0 / (1.987 * T_K)) * pow(sg, 5.0);
    const double tem7 = 9.53e21  * exp(-282700.0 / (1.987 * T_K))
                      * pow(sg, 10.0) * pow(time_s, 3.0);
    const double TP = (tem1 + tem2 + tem3) * tem4;
    const double TS = tem5 + tem6;
    const double TT = 4.0 * tem7;
    const double total = (ircreep + TP + TS + TT) * 3.6e3 * 0.01; // %/h → /s scaling
    return min(total, 1.0e-4) / 3.6e3; // cap and convert /h → /s
}

// Void swelling volumetric strain [-]  — Cswel.for model 1, Hofmann (1985)
//   neuflue : neutron fluence [10^22 n/cm^2]
double VoidSwelling(double neuflue, double T_K) {
    const double ttt  = T_K - 273.15;
    const double R    = 0.085 * exp(-1.0e-4 * (ttt - 400.0) * (ttt - 400.0));
    const double D    = 0.01 * 0.15 * (1.0 - exp(-0.1 * neuflue));
    const double tem1 = (1.0 + exp(0.75 * (14.2 - neuflue)))
                       / (1.0 + exp(0.75 * 14.2));
    const double s0   = 0.01 * R * (neuflue + 1.0 / 0.75 * log(tem1));
    return s0 + D;
}

// Burst stress [Pa]  — Csigb.for (HT-9 branch, bug-corrected)
// Fortran: ht9_tem = (tk-tk-200)/200  where tk-tk = 0 → always -1 (bug).
// Corrected: ht9_tem = (T_Celsius - 200) / 200
double BurstStress(double T_K, double ssy0) {
    const double T_C    = T_K - 273.15;
    const double factor = std::tanh((T_C - 200.0) / 200.0);
    return ssy0 * (1.1 - 0.1 * factor);
}

// Void swelling RATE [1/s] — Cswel.for icswel=2, SAS4A User's Manual
// piecewise-linear-in-T coefficient table.  The dose >= 100 dpa onset gate
// is applied generically by the caller (FredMNaSolver); this function
// returns the material's own T-dependent rate curve unconditionally.
double SAS4ASwellingRate(double T_K, double flux_1e22, double dconv) {
    double coeff;
    if (T_K <= 673.0)
        coeff = 8.33e-5 + (1.0e-4  - 8.33e-5) * (T_K - 623.0) / 50.0;
    else if (T_K <= 723.0)
        coeff = 1.0e-4  + (8.33e-5 - 1.0e-4)  * (T_K - 673.0) / 50.0;
    else if (T_K <= 773.0)
        coeff = 8.33e-5 + (5.0e-5  - 8.33e-5) * (T_K - 723.0) / 50.0;
    else if (T_K <= 823.0)
        coeff = 7.5e-5  + (2.5e-5  - 7.5e-5)  * (T_K - 773.0) / 50.0;
    else if (T_K <= 873.0)
        coeff = 5.0e-5  + (1.25e-5 - 5.0e-5)  * (T_K - 823.0) / 50.0;
    else
        return 0.0;
    return coeff * flux_1e22 * dconv;
}

// Cladding wastage correlations — Clanth.for (precipitation kinetics) and
// MFUEL App. 9.3 (lanthanide tracking).  Diffusion coefficients are pure
// Arrhenius; the shared irradiation-enhanced-diffusivity floor (T_floor,
// ~800 K) is applied generically by the scheme code in
// FredMNaCladWastage.hpp before calling these, not here.
constexpr double kWastagePkD0   = 1632573.01020072;  // [m2/s] Clanth.for fit
constexpr double kWastagePkQ    = 333962.40004088;   // [J/mol]
constexpr double kWastageLaD0   = 5.0;               // [m2/s] MFUEL App. 9.3
constexpr double kWastageLaQ    = 250.0e3;           // [J/mol]
constexpr double kWastageUSol   = 6.4e27;            // [#/m3] HT-9 (D9 would be 0.5e27)

double WastageDiffusionPK(double T_K) {
    return kWastagePkD0 * exp(-kWastagePkQ / (8.314 * T_K));
}
double WastageDiffusionLaTracking(double T_K) {
    return kWastageLaD0 * exp(-kWastageLaQ / (8.314 * T_K));
}

} // anonymous namespace

// ---------------------------------------------------------------------------
namespace fred {

HT9::HT9(double rho0) : m_rho0(rho0) {}

double HT9::thermalConductivity(double T) const {
    return ThermalConductivity(T);
}

double HT9::heatCapacity(double T) const {
    return HeatCapacity(T);
}

double HT9::thermalExpansionStrain(double T) const {
    return ThermalExpansionStrain(T);
}

// Returns E in MPa (CladdingMaterial interface unit).
double HT9::youngsModulus(double T) const {
    return ElasticModulus(T) * 1.0e-6;
}

// Returns mid-temperature average; use the static PoissonRatio helper
// directly when temperature-dependence is needed.
double HT9::poissonRatio() const {
    return 0.28;
}

double HT9::meyerHardness(double T) const {
    return MeyerHardness(T);
}

double HT9::referenceDensity() const { return m_rho0; }

// HT9-specific extended methods (beyond CladdingMaterial interface).

double HT9::yieldStress(double T) const {
    return YieldStress(T);
}

// sig_Pa  : effective stress [Pa]
// qqv     : volumetric power density [W/m3] (neutron flux proxy)
// time_s  : elapsed time [s]
double HT9::creepRate(double T, double sig_Pa, double qqv, double time_s) const {
    return CreepRate(T, sig_Pa, qqv, time_s);
}

// neuflue : neutron fluence [10^22 n/cm^2]
double HT9::voidSwelling(double neuflue, double T) const {
    return VoidSwelling(neuflue, T);
}

// ssy0 : yield stress [Pa] at the same temperature (from yieldStress(T))
double HT9::burstStress(double T_K, double ssy0) const {
    return BurstStress(T_K, ssy0);
}

double HT9::voidSwellingSAS4ARate(double T_K, double flux_1e22,
                                   double dconv, double dose_dpa) const {
    (void)dose_dpa;  // onset gate applied by the caller
    return SAS4ASwellingRate(T_K, flux_1e22, dconv);
}

double HT9::wastageDiffusionPK(double T_K) const {
    return WastageDiffusionPK(T_K);
}
double HT9::wastageSolubilityFuel() const { return 0.0; }
double HT9::wastageSolubilityClad() const { return 0.1; }
double HT9::wastageDiffusionLaTracking(double T_K) const {
    return WastageDiffusionLaTracking(T_K);
}
double HT9::wastageUptakeCapacity() const { return kWastageUSol; }

} // namespace fred
