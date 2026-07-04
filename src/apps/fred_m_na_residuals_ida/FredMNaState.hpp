#pragma once
#include "platform/FuelRodState.hpp"
#include "apps/fred_m_na/FredMNaGrsis.hpp"
#include <vector>
#include <string>

namespace fred {

// Per-fuel-node irradiation state for FRED-M-Na.
struct FredMNaNodeState {
    // Composition (updated by Zr redistribution)
    double zr_wf = 0.0;   // Zr weight fraction [-]
    double pu_wf = 0.0;   // Pu weight fraction [-]
    double ur_wf = 0.0;   // U weight fraction [-]
    double zr_at = 0.0;   // Zr atom fraction [-]
    double pu_at = 0.0;   // Pu atom fraction [-]
    double ur_at = 0.0;   // U atom fraction [-]
    double c_zr  = 0.0;   // Zr atomic density [atom/m3]
    double mass  = 0.0;   // node mass [kg] (fixed after init)
    double dvol  = 0.0;   // node volume [m3] (fixed after init)

    // Phase and sodium infiltration
    std::string phase = "alpha"; // current phase
    double pfrac = 1.0;          // phase fraction
    double psod  = 0.0;          // sodium infiltration fraction [-]

    // Porosity state (needed for thermal conductivity)
    double poros_tot = 0.0;  // total porosity [-]  (driven by swelling)
    double poros_gas = 0.0;  // gas-filled porosity [-]

    // GRSIS bubble-swelling state (updated in FredMNaSolver::afterAcceptedStep)
    GrsisFuelNodeState grsis;
};

// Per-axial-layer state for FRED-M-Na.
//
// Extends AxialLayerState with irradiation fields:
//   bup        : local burnup [MWd/kgU]
//   bup_FIMA   : local burnup [FIMA] (from Daniele formula)
//   buhard_FIMA: burnup at hard-contact event [FIMA]
//   flag       : gap contact state: "open", "soft", or "clos"
//   efs[nf]    : volumetric fuel swelling strain [-]
//   nodes[nf]  : per-node irradiation state
//   xwast      : cladding wastage thickness [m]
//   clfuel     : lanthanide concentration in fuel [-]
//   ntot       : fuel atom density [atom/m3]
struct FredMNaLayerState : public AxialLayerState {
    double bup        = 0.0;
    double dbup       = 0.0;
    double bup_FIMA   = 0.0;
    double buhard_FIMA= 0.0;

    std::string flag  = "open"; // "open" / "soft" / "clos"

    // NOTE: volumetric swelling strain efs[] is inherited from AxialLayerState.
    //       Do not redeclare it here — StressStrain reads AxialLayerState::efs.
    std::vector<double> defs;  // swelling rate [1/s], size nf
    std::vector<double> dec;   // cladding hoop creep rate [1/s], size nc

    // d(efs)/dt implied by GRSIS's d(swtot)/dt / 3 (volumetric -> linear strain),
    // updated once per accepted step in FredMNaSolver::afterAcceptedStep. Used by
    // computeIrradiationResiduals in place of the solid-FP-only empirical rate
    // once GRSIS has produced a value (nodes[i].grsis.swtot > 0), so gap closure
    // is driven by GRSIS's total (solid + gas bubble) swelling, not just solid FPs.
    std::vector<double> defs_grsis;  // [1/s], size nf

    // Directional swelling strain accumulators for Baseir.for anisotropy:
    // soft  → Δefsz=0, Δefsh=Δefsr=0.4995·ΔSwtot
    // open/clos → Δefsz=Δefsh=Δefsr=max(0, ΔSwtot/3)
    // Consumed by FredMNaStressStrain; NOT in AxialLayerState (platform-untouched).
    std::vector<double> efsz;  // directional swelling: axial  [nf]
    std::vector<double> efsh;  // directional swelling: hoop   [nf]
    std::vector<double> efsr;  // directional swelling: radial [nf]

    std::vector<FredMNaNodeState> nodes; // size nf

    // Cladding wastage (tracked per layer)
    double xwast  = 0.0;  // wastage thickness [m]
    double clfuel = 0.0;  // lanthanide concentration [-]
    double ntot   = 0.0;  // fuel atom density [atom/m3]

    explicit FredMNaLayerState(int nf_, int nc_)
        : AxialLayerState(nf_, nc_),
          defs(nf_, 0.0),
          dec (nc_, 0.0),
          defs_grsis(nf_, 0.0),
          efsz(nf_, 0.0),
          efsh(nf_, 0.0),
          efsr(nf_, 0.0),
          nodes(nf_)
    {}
};

// Full rod state for FRED-M-Na.
struct FredMNaRodState {
    int nf, nc, nz;

    // Global gas state
    double gpres       = 0.1;  // internal gas pressure [MPa]
    double mu0         = 0.0;  // initial fill-gas inventory [mol]
    double fggen       = 0.0;  // total fission gas generated [mol]
    double fgrel       = 0.0;  // total fission gas released [mol]
    double dfggen      = 0.0;  // fission gas generation rate [mol/s]
    double dfgrel      = 0.0;  // fission gas release rate [mol/s]
    double bupave_FIMA = 0.0;  // spatially averaged burnup [FIMA]

    // Fuel axial anisotropy factor (rod-scalar, see upuzrFanis / Baseir.for).
    // Drives the "soft" contact criterion: gap - fanis_coef*gap0 <= rough.
    // Frozen once the peak-power layer's burnup exceeds 0.5% FIMA.
    double fanis_F    = 0.0;
    double fanis_coef = 0.0;

    std::vector<FredMNaLayerState> layers;

    FredMNaRodState(int nf_, int nc_, int nz_)
        : nf(nf_), nc(nc_), nz(nz_),
          layers(nz_, FredMNaLayerState(nf_, nc_))
    {}
};

} // namespace fred
