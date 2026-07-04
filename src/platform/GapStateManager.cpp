#include "GapStateManager.hpp"
#include <iostream>

namespace fred {

GapStateManager::GapStateManager(const FuelRodGeometry& geom)
{
    const bool fab_open = (geom.rci0 - geom.rfo0) > (geom.ruff + geom.rufc);
    m_gapOpen.assign(geom.nz, fab_open);
    m_axialOffset.assign(geom.nz, 0.0);
}

void GapStateManager::setGapOpen(int layer, bool open) {
    m_gapOpen[layer] = open;
    if (open) m_axialOffset[layer] = 0.0;
}

double GapStateManager::applyGapClosed(int layer, double efz, double ez) {
    const double offset = efz - ez;
    m_axialOffset[layer] = offset;
    m_gapOpen[layer]     = false;
    std::cout << "  Axial offset at closure (layer " << layer << "): " << offset << "\n";
    return offset;
}

void GapStateManager::applyToLayer(int layer, AxialLayerState& s) const {
    s.gapOpen     = m_gapOpen[layer];
    s.axialOffset = m_axialOffset[layer];
}

} // namespace fred
