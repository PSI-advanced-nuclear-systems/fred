#pragma once
#include "platform/SubchannelMode.hpp"
#include "platform/SubchannelCoolantProperties.hpp"

namespace fred {

// Subchannel coolant field updater for FRED-M-Na (sodium fast reactor).
// Implements axial energy balance + Peclet-number HTC correlations.
// Reference: Baseir.for (FRED-M Timpano), lines ~270-287.
class FredMNaSubchannelMode : public SubchannelMode {
public:
    enum class HtcCorrelation { kMikityuk, kSubbotin };

    FredMNaSubchannelMode(double dhyd, double xarea, double flowr,
                           TimeSeries T_inlet_fn,
                           double rfo0,
                           const SubchannelCoolantProperties& props,
                           int nz,
                           HtcCorrelation corr = HtcCorrelation::kMikityuk);

    void updateCoolantField(double t,
                            const double* qqv,
                            const double* dz) override;
private:
    double m_rfo0;
    const SubchannelCoolantProperties& m_props;
    HtcCorrelation m_corr;
};

} // namespace fred
