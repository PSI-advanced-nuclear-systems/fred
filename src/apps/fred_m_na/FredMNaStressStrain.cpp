#include "FredMNaStressStrain.hpp"
#include "platform/Constants.hpp"
#include "platform/FuelPelletMaterial.hpp"
#include "platform/CladdingMaterial.hpp"
#include <cmath>
#include <cassert>

namespace fred {

FredMNaStressStrain::FredMNaStressStrain(const FuelRodGeometry& geom,
                                          const FuelPelletMaterial& fuel,
                                          const CladdingMaterial& clad)
    : m_geom(geom), m_fuel(fuel), m_clad(clad)
{}

int FredMNaStressStrain::neq() const {
    // Same layout/count as StressStrain::neq() — depends only on nf/nc.
    return 6 * m_geom.nf + 6 * m_geom.nc + 2;
}

void FredMNaStressStrain::computeResiduals(const FredMNaLayerState& s,
                                            double gpres, double pcool,
                                            std::vector<double>& r) const
{
    const int nf = m_geom.nf;
    const int nc = m_geom.nc;
    assert((int)r.size() >= neq());

    const bool open_ = (s.flag == "open");
    const bool soft_  = (s.flag == "soft");

    // Material moduli at current temperatures
    std::vector<double> Ef(nf), nuF(nf);
    for (int i = 0; i < nf; ++i) {
        Ef[i]  = m_fuel.youngsModulus(s.T[i], m_fuel.referenceDensity());
        nuF[i] = m_fuel.poissonRatio();
    }
    std::vector<double> Ec(nc), nuC(nc);
    for (int i = 0; i < nc; ++i) {
        Ec[i]  = m_clad.youngsModulus(s.T[nf + i]);
        nuC[i] = m_clad.poissonRatio();
    }

    // Average strains at ring boundaries (used in compatibility and equilibrium)
    // efr_half[i] = 0.5*(efr[i] + efr[i+1])  for fuel rings
    std::vector<double> efr_half(nf-1), sigfh_half(nf-1);
    for (int i = 0; i < nf-1; ++i) {
        efr_half[i]  = 0.5 * (s.efr[i] + s.efr[i+1]);
        sigfh_half[i]= 0.5 * (s.sigfh[i] + s.sigfh[i+1]);
    }
    std::vector<double> er_half(nc-1), sigh_half(nc-1);
    for (int i = 0; i < nc-1; ++i) {
        er_half[i]  = 0.5 * (s.er[i] + s.er[i+1]);
        sigh_half[i]= 0.5 * (s.sigh[i] + s.sigh[i+1]);
    }

    int n = 0; // local equation counter

    // =========================================================================
    // FUEL MECHANICS
    // =========================================================================

    // --- Fuel thermal expansion AE ---
    for (int i = 0; i < nf; ++i) {
        double eft_target = m_fuel.thermalExpansionStrain(s.T[i]);
        r[n++] = s.eft[i] - eft_target;
    }

    // --- Fuel Hooke's law (thermo-elastic + directional swelling) ---
    // Elastic strain = total strain - thermal expansion - directional swelling.
    // efsz/efsh/efsr are the Baseir.for anisotropic swelling accumulators
    // (isotropic ΔSwtot/3 in open/clos, axial-suppressed 0.4995*ΔSwtot in soft
    // for hoop/radial only), updated per-step by FredMNaGapBehavior — this
    // replaces StressStrain's hardcoded s.efs[i]/3.0 isotropic term.
    // Hoop:   E * (efh - eft - efsh) = sigfh - nu*(sigfr + sigfz)
    for (int i = 0; i < nf; ++i) {
        double eel = s.efh[i] - s.eft[i] - s.efsh[i];
        r[n++] = Ef[i]*eel - s.sigfh[i] + nuF[i]*(s.sigfr[i] + s.sigfz[i]);
    }
    // Radial: E * (efr - eft - efsr) = sigfr - nu*(sigfh + sigfz)
    for (int i = 0; i < nf; ++i) {
        double eel = s.efr[i] - s.eft[i] - s.efsr[i];
        r[n++] = Ef[i]*eel - s.sigfr[i] + nuF[i]*(s.sigfh[i] + s.sigfz[i]);
    }
    // Axial:  E * (efz - eft - efsz) = sigfz - nu*(sigfh + sigfr)
    for (int i = 0; i < nf; ++i) {
        double eel = s.efz - s.eft[i] - s.efsz[i];   // efz is uniform over the layer
        r[n++] = Ef[i]*eel - s.sigfz[i] + nuF[i]*(s.sigfh[i] + s.sigfr[i]);
    }

    // --- Fuel strain compatibility: (rr-rl)*efr_avg + rl*efh[l] - rr*efh[r] = 0 ---
    for (int i = 0; i < nf-1; ++i) {
        r[n++] = s.drf[i] * efr_half[i]
                 + s.rad[i]   * s.efh[i]
                 - s.rad[i+1] * s.efh[i+1];
    }

    // --- Fuel stress equilibrium: (rr-rl)*sigfh_avg + rl*sigfr[l] - rr*sigfr[r] = 0 ---
    for (int i = 0; i < nf-1; ++i) {
        r[n++] = s.drf[i] * sigfh_half[i]
                 + s.rad[i]   * s.sigfr[i]
                 - s.rad[i+1] * s.sigfr[i+1];
    }

    // --- Fuel BC-1: global axial force balance ---
    // Legacy: soft behaves like clos (shared axial load) — item 2 of the
    // Baseir.for physics table.
    {
        double sum = 0.0;
        for (int i = 0; i < nf; ++i)
            sum += m_geom.area0[i] * s.sigfz[i];

        if (open_) {
            // Gap open: fuel axial stress balanced by gas pressure only
            sum += m_geom.azf0 * gpres;
        } else {
            // soft or clos: combined fuel+clad axial equilibrium.
            for (int i = 0; i < nc; ++i)
                sum += m_geom.area0[nf + i] * s.sigz[i];
            sum += pcool * PI * m_geom.rco0 * m_geom.rco0;
        }
        r[n++] = sum;
    }

    // --- Fuel BC-2: inner radial boundary (unrelated to gap-contact state) ---
    if (m_geom.rfi0 == 0.0) {
        r[n++] = s.sigfr[0] - s.sigfh[0];
    } else {
        r[n++] = s.sigfr[0] + gpres;
    }

    // --- Fuel BC-3: outer radial / gap interface ---
    // Legacy: soft behaves like clos (strain continuity) — item 3 of the
    // Baseir.for physics table.
    if (open_) {
        r[n++] = s.sigfr[nf-1] + gpres;
    } else {
        // soft or clos: preserve axial strain offset frozen at closure.
        r[n++] = s.efz - s.ez - s.axialOffset;
    }

    // =========================================================================
    // CLADDING MECHANICS
    // =========================================================================

    // --- Clad thermal expansion AE ---
    for (int i = 0; i < nc; ++i) {
        double et_target = m_clad.thermalExpansionStrain(s.T[nf + i]);
        r[n++] = s.et[i] - et_target;
    }

    // --- Clad Hooke's law (thermo-elastic + creep) ---
    for (int i = 0; i < nc; ++i) {
        double eel = s.eh[i] - s.et[i] - s.ec[i];
        r[n++] = Ec[i]*eel - s.sigh[i] + nuC[i]*(s.sigr[i] + s.sigz[i]);
    }
    for (int i = 0; i < nc; ++i) {
        double eel = s.er[i] - s.et[i] + 0.5*s.ec[i];
        r[n++] = Ec[i]*eel - s.sigr[i] + nuC[i]*(s.sigh[i] + s.sigz[i]);
    }
    for (int i = 0; i < nc; ++i) {
        double eel = s.ez - s.et[i] + 0.5*s.ec[i];
        r[n++] = Ec[i]*eel - s.sigz[i] + nuC[i]*(s.sigh[i] + s.sigr[i]);
    }

    // --- Clad strain compatibility ---
    for (int i = 0; i < nc-1; ++i) {
        r[n++] = s.drc[i] * er_half[i]
                 + s.rad[nf+i]   * s.eh[i]
                 - s.rad[nf+i+1] * s.eh[i+1];
    }

    // --- Clad stress equilibrium ---
    for (int i = 0; i < nc-1; ++i) {
        r[n++] = s.drc[i] * sigh_half[i]
                 + s.rad[nf+i]   * s.sigr[i]
                 - s.rad[nf+i+1] * s.sigr[i+1];
    }

    // --- Clad BC-1: axial or interface ---
    // Legacy: soft behaves like clos (cladding axial BC / interface
    // continuity) — item 4 of the Baseir.for physics table.
    if (open_) {
        double sum = 0.0;
        for (int i = 0; i < nc; ++i)
            sum += m_geom.area0[nf + i] * s.sigz[i];
        sum += pcool * PI * m_geom.rco0 * m_geom.rco0
               - gpres * PI * m_geom.rci0 * m_geom.rci0;
        r[n++] = sum;
    } else {
        // soft or clos: radial stress continuity at fuel–clad interface.
        r[n++] = s.sigfr[nf-1] - s.sigr[0];
    }

    // --- Clad BC-2: inner radial surface ---
    // Legacy: soft behaves like open (gas still surrounds the cladding ID
    // because contact is localized, not full annular) — item 5 of the
    // Baseir.for physics table. This is the one BC site where soft groups
    // with open rather than clos.
    if (open_ || soft_) {
        r[n++] = s.sigr[0] + gpres;
    } else {
        // clos: geometric contact constraint (gap = roughness)
        r[n++] = s.gap - (m_geom.ruff + m_geom.rufc);
    }

    // --- Clad BC-3: outer radial surface (coolant pressure) ---
    r[n++] = s.sigr[nc-1] + pcool;

    assert(n == neq());
}

} // namespace fred
