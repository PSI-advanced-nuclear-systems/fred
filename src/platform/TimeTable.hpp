#pragma once
#include <vector>

namespace fred {

// -----------------------------------------------------------------------
// TimeTable — piecewise-linear interpolation of a scalar vs time.
//
// operator()(t) clamps to the first/last value outside [t_min, t_max].
// Shared by all FRED solver drivers (ROD, OX, M-Na).
// -----------------------------------------------------------------------
class TimeTable {
public:
    TimeTable() = default;
    TimeTable(std::vector<double> times, std::vector<double> values);

    double operator()(double t) const;
    bool   empty() const { return m_t.empty(); }

    // Public for lambda captures inside solver drivers.
    std::vector<double> m_t, m_v;
};

} // namespace fred
