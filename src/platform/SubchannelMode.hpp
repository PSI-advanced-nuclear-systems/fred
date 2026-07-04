#pragma once
#include <functional>
#include <vector>

namespace fred {

// Piecewise-linear time-series type for boundary conditions.
using TimeSeries = std::function<double(double)>;

// Abstract base for subchannel coolant field update strategies.
// Concrete subclasses (e.g. FredMNaSubchannelMode) implement updateCoolantField
// for a specific coolant type and HTC correlation.
class SubchannelMode {
public:
    SubchannelMode(double dhyd, double xarea, double flowr,
                   TimeSeries T_inlet_fn, int nz);
    virtual ~SubchannelMode() = default;

    // Called by RodResiduals before the per-layer residual loop.
    // qqv[nz]: volumetric power per layer [W/m3]. dz[nz]: layer heights [m].
    virtual void updateCoolantField(double t,
                                    const double* qqv,
                                    const double* dz) = 0;

    double T_co(int j) const { return m_T_co[j]; }
    double htc (int j) const { return m_htc[j];  }

protected:
    double     m_dhyd, m_xarea, m_flowr;
    TimeSeries m_T_inlet_fn;
    std::vector<double> m_T_co;
    std::vector<double> m_htc;
};

} // namespace fred
