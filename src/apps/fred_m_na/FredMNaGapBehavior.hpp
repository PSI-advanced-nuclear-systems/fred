#pragma once
#include "apps/fred_m_na/FredMNaState.hpp"

namespace fred {

// Owns the gap-contact ratchet and the per-step directional swelling update
// for FRED-M-Na (U-Pu-Zr metallic fuel).
//
// Ratchet (monotonic, legacy Baseir.for, never reopens for upuzr):
//   open → soft  when  gap - fanis_coef * gap0  <=  ruff + rufc
//   soft → clos  when  gap  <=  ruff + rufc  (gapOpen == false)
//
// Directional swelling (Baseir.for lines 606-627):
//   soft:      Δefsz = 0,             Δefsh = Δefsr = 0.4995 * ΔSwtot
//   open/clos: Δefsz = Δefsh = Δefsr = max(0, ΔSwtot / 3)
//
// Called from FredMNaSolver::afterAcceptedStep once per accepted fixed
// (backward-Euler) step -- FredMNaSolver does not use SUNDIALS IDA.
// \coderef{fred_m_na_gap_behavior}
class FredMNaGapBehavior {
public:
    // Update s.flag (ratchet) and s.efsz/efsh/efsr (swelling partition) for
    // one axial layer.
    //
    //   s           : layer state to update (flag, efsz, efsh, efsr)
    //   gap0        : as-fabricated gap [m]  (rci0 - rfo0)
    //   rough       : combined surface roughness [m]  (ruff + rufc)
    //   fanis_coef  : soft-contact trigger coefficient (from upuzrFanis)
    //   swtot_new[] : current total volumetric swelling per fuel node [nf]
    //   swtot_old[] : total volumetric swelling at previous step    [nf]
    static void update(FredMNaLayerState& s,
                       double gap0, double rough, double fanis_coef,
                       const double* swtot_new, const double* swtot_old,
                       int nf);
};

} // namespace fred
