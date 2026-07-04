#include "TimeTable.hpp"
#include <algorithm>

namespace fred {

TimeTable::TimeTable(std::vector<double> times, std::vector<double> values)
    : m_t(std::move(times)), m_v(std::move(values))
{}

double TimeTable::operator()(double t) const {
    if (m_t.empty())         return 0.0;
    if (t <= m_t.front())    return m_v.front();
    if (t >= m_t.back())     return m_v.back();
    auto it  = std::lower_bound(m_t.begin(), m_t.end(), t);
    int  idx = (int)(it - m_t.begin()) - 1;
    double f = (t - m_t[idx]) / (m_t[idx+1] - m_t[idx]);
    return m_v[idx] + f * (m_v[idx+1] - m_v[idx]);
}

} // namespace fred
