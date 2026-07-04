#include "apps/fred_m_na/FredMNaGapBehavior.hpp"
#include <algorithm>

namespace fred {

void FredMNaGapBehavior::update(FredMNaLayerState& s,
                                 double gap0, double rough, double fanis_coef,
                                 const double* swtot_new, const double* swtot_old,
                                 int nf)
{
    // --- Gap-contact ratchet (monotonic: never reopens for upuzr) ---
    if (s.flag != "clos") {
        if (!s.gapOpen) {
            s.flag = "clos";
        } else if (s.gap - fanis_coef * gap0 <= rough) {
            s.flag = "soft";
        }
        // else: stays "open"
    }

    // --- Directional swelling update (Baseir.for lines 606-627) ---
    for (int i = 0; i < nf; ++i) {
        const double dsw = swtot_new[i] - swtot_old[i];
        if (s.flag == "soft") {
            s.efsz[i] += 0.0;
            s.efsh[i] += 0.4995 * dsw;
            s.efsr[i] += 0.4995 * dsw;
        } else {
            const double iso = std::max(0.0, dsw / 3.0);
            s.efsz[i] += iso;
            s.efsh[i] += iso;
            s.efsr[i] += iso;
        }
    }
}

bool FredMNaGapBehavior::stateChanged(const std::string& before,
                                       const std::string& after)
{
    return before != after;
}

} // namespace fred
