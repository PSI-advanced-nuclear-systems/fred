#include "SubchannelMode.hpp"

namespace fred {

SubchannelMode::SubchannelMode(double dhyd, double xarea, double flowr,
                                TimeSeries T_inlet_fn, int nz)
    : m_dhyd(dhyd), m_xarea(xarea), m_flowr(flowr),
      m_T_inlet_fn(std::move(T_inlet_fn)),
      m_T_co(nz, m_T_inlet_fn(0.0)),
      m_htc(nz, 0.0)
{}

} // namespace fred
