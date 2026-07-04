#include "FredMNaSubchannelMode.hpp"
#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace fred {

FredMNaSubchannelMode::FredMNaSubchannelMode(
        double dhyd, double xarea, double flowr,
        TimeSeries T_inlet_fn, double rfo0,
        const SubchannelCoolantProperties& props,
        int nz, HtcCorrelation corr)
    : SubchannelMode(dhyd, xarea, flowr, std::move(T_inlet_fn), nz),
      m_rfo0(rfo0), m_props(props), m_corr(corr)
{}

void FredMNaSubchannelMode::updateCoolantField(double t,
                                                const double* qqv,
                                                const double* dz)
{
    const int nz = static_cast<int>(m_T_co.size());
    const double T_inlet = m_T_inlet_fn(t);

    // Pitch-to-diameter ratio fixed to ESFR-SMART assembly geometry (Baseir.for xx=1.1)
    const double xx = 1.1;

    // Layer rod powers [W]: P = qqv * pi * rfo^2 * dz
    std::vector<double> P(nz);
    for (int j = 0; j < nz; ++j)
        P[j] = qqv[j] * M_PI * m_rfo0 * m_rfo0 * dz[j];

    for (int j = 0; j < nz; ++j) {
        const double T_prev = (j == 0) ? T_inlet : m_T_co[j - 1];
        const double Cp  = m_props.cp(T_prev);
        const double rho = m_props.rho(T_prev);
        const double k   = m_props.k(T_prev);

        // Axial energy balance (Baseir.for lines 215-219)
        if (j == 0)
            m_T_co[0] = T_inlet + 0.5 * P[0] / (m_flowr * Cp);
        else
            m_T_co[j] = m_T_co[j-1] + 0.5 * (P[j-1] + P[j]) / (m_flowr * Cp);

        // Coolant velocity and Peclet number
        const double vcool = m_flowr / (rho * m_xarea);
        const double Pe    = vcool * m_dhyd * rho * Cp / k;

        // HTC correlation (Baseir.for lines 281-287)
        double htc_val;
        if (m_corr == HtcCorrelation::kMikityuk) {
            // Mikityuk: liquid-metal rod-bundle correlation
            htc_val = k / m_dhyd * (0.047 * (1.0 - std::exp(-3.8 * (xx - 1.0)))
                                     * (std::pow(Pe, 0.77) + 250.0));
        } else {
            // Subbotin: simplified liquid-metal correlation
            htc_val = k / m_dhyd * (5.0 + 0.025 * std::pow(Pe, 0.8));
        }
        m_htc[j] = htc_val;
    }
}

} // namespace fred
