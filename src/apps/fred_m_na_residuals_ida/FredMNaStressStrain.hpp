#pragma once
#include "platform/FuelRodGeometry.hpp"
#include "apps/fred_m_na/FredMNaState.hpp"
#include <vector>

namespace fred {

class FuelPelletMaterial;
class CladdingMaterial;

// FredMNaStressStrain — FRED-M-Na's own three-state (open/soft/clos) mechanics
// residual class (Decision D2/D8 of the gap-behaviour-model refactor).
//
// This is a .cpp-only fork of platform/StressStrain.cpp: StressStrain's
// computeResiduals() is non-virtual and its private internals are not
// reachable via composition, so this class copies the full residual layout
// and modifies exactly two things relative to the original:
//
//   1. The 4 BC branch sites that read AxialLayerState::gapOpen (a bool) now
//      read FredMNaLayerState::flag (a 3-way "open"/"soft"/"clos" string),
//      per the legacy Baseir.for physics table:
//        Fuel BC-1 (axial force balance)   : open vs {soft,clos}
//        Fuel BC-3 (outer radial)          : open vs {soft,clos}
//        Clad BC-1 (axial/interface)       : open vs {soft,clos}
//        Clad BC-2 (inner radial)          : {open,soft} vs clos
//   2. The isotropic swelling term efs[i]/3.0 is replaced by the directional
//      accumulators efsz[i]/efsh[i]/efsr[i] (Baseir.for anisotropy partition,
//      updated per-step by FredMNaGapBehavior).
//
// platform/StressStrain.{hpp,cpp} are completely untouched (Decision D1);
// FRED-ROD/FRED-OX output cannot be affected by this class.
// \coderef{fred_m_na_stress_strain}
class FredMNaStressStrain {
public:
    FredMNaStressStrain(const FuelRodGeometry& geom,
                         const FuelPelletMaterial& fuel,
                         const CladdingMaterial& clad);

    // Same equation count/layout as StressStrain::neq() (depends only on nf/nc).
    int neq() const;

    // Fill residual vector r[0..neq()-1] for axial layer with state `s`.
    // \coderef{fred_m_na_stress_strain_residual}
    void computeResiduals(const FredMNaLayerState& s,
                          double gpres, double pcool,
                          std::vector<double>& r) const;

private:
    const FuelRodGeometry&    m_geom;
    const FuelPelletMaterial& m_fuel;
    const CladdingMaterial&   m_clad;
};

} // namespace fred
