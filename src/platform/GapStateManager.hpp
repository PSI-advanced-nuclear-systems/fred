#pragma once
#include "FuelRodGeometry.hpp"
#include "FuelRodState.hpp"
#include <vector>

namespace fred {

// Tracks per-layer gap open/closed state and the axial strain offset frozen at
// closure. Shared across all residual assemblers to avoid duplicating this logic.
//
// Closure offset: when the gap closes, C = efz - ez is frozen so the
// constraint  efz - ez = C  (equivalent to defz = dez) is enforced at
// subsequent IDA steps. The offset is cleared when the gap reopens.
// \coderef{gap_state_manager}
class GapStateManager {
public:
    // Initialise from as-fabricated geometry. Sets gapOpen = (fab gap > roughness sum).
    GapStateManager(const FuelRodGeometry& geom);

    // Called by the solver on root events.
    void setGapOpen(int layer, bool open);

    // Called at gap CLOSURE with efz and ez unpacked from the current IDA state.
    // Records axialOffset = efz - ez and returns the stored offset.
    double applyGapClosed(int layer, double efz, double ez);

    // Propagate current gap state into a layer state struct.
    void applyToLayer(int layer, AxialLayerState& s) const;

    bool   isGapOpen(int layer)   const { return m_gapOpen[layer]; }
    double axialOffset(int layer) const { return m_axialOffset[layer]; }
    void setAxialOffset(int layer, double offset) { m_axialOffset[layer] = offset; }

private:
    std::vector<bool>   m_gapOpen;
    std::vector<double> m_axialOffset;
};

} // namespace fred
