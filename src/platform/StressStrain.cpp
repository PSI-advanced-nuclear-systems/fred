#include "StressStrain.hpp"
#include "Constants.hpp"
#include "FuelPelletMaterial.hpp"
#include "CladdingMaterial.hpp"
#include <cmath>
#include <cassert>

namespace fred {

StressStrain::StressStrain(const FuelRodGeometry& geom,
                           const FuelPelletMaterial& fuel,
                           const CladdingMaterial& clad)
    : m_geom(geom), m_fuel(fuel), m_clad(clad)
{}

int StressStrain::neq() const {
    // 4*nf + 4*nc: thermal expansion AE + Hooke's law (theta, r, z) per node
    // + 2*(nf-1) + 2*(nc-1): strain compatibility + stress equilibrium per ring pair
    // 6 boundary conditions (emitted in this order):
    // Fuel BC-1: open: fuel axial force balance;  closed: combined fuel+clad axial balance
    // Fuel BC-2: solid: sigfr[0]=sigfh[0];  annular: sigfr[0]=-gpres
    // Fuel BC-3: open: sigfr[nf-1]=-gpres;  closed: efz_fuel=efz_clad
    // Clad BC-1: open: clad axial force balance;  closed: traction continuity sigfr[nf-1]=sigr[0]
    // Clad BC-2: open: sigr[0]=-gpres;  closed: gap=roughness
    // Clad BC-3: sigr[nc-1]=-pcool
    return 6 * m_geom.nf + 6 * m_geom.nc + 2;
}

void StressStrain::computeResiduals(const AxialLayerState& s,
                                    double gpres, double pcool,
                                    std::vector<double>& r) const
{
    const int nf = m_geom.nf;
    const int nc = m_geom.nc;
    assert((int)r.size() >= neq());

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
    // eft[i] = ftexp(T[i]) - ftexp(T_REF)
    for (int i = 0; i < nf; ++i) {
        double eft_target = m_fuel.thermalExpansionStrain(s.T[i]);
        r[n++] = s.eft[i] - eft_target;
    }

    // --- Fuel Hooke's law (thermo-elastic + swelling) ---
    // Elastic strain = total strain - thermal expansion - isotropic swelling (efs/3 per direction)
    // efs is the volumetric swelling strain (0 for FRED-ROD, irradiation swelling for FRED-OX).
    // Hoop:   E * (efh - eft - efs/3) = sigfh - nu*(sigfr + sigfz)
    for (int i = 0; i < nf; ++i) {
        double eel = s.efh[i] - s.eft[i] - s.efs[i] / 3.0;
        r[n++] = Ef[i]*eel - s.sigfh[i] + nuF[i]*(s.sigfr[i] + s.sigfz[i]);
    }
    // Radial: E * (efr - eft - efs/3) = sigfr - nu*(sigfh + sigfz)
    for (int i = 0; i < nf; ++i) {
        double eel = s.efr[i] - s.eft[i] - s.efs[i] / 3.0;
        r[n++] = Ef[i]*eel - s.sigfr[i] + nuF[i]*(s.sigfh[i] + s.sigfz[i]);
    }
    // Axial:  E * (efz - eft - efs/3) = sigfz - nu*(sigfh + sigfr)
    for (int i = 0; i < nf; ++i) {
        double eel = s.efz - s.eft[i] - s.efs[i] / 3.0;   // efz is uniform over the layer
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
    {
        double sum = 0.0;
        for (int i = 0; i < nf; ++i)
            sum += m_geom.area0[i] * s.sigfz[i];

        if (s.gapOpen) {
            // Gap open: fuel axial stress balanced by gas pressure only
            sum += m_geom.azf0 * gpres;
        } else {
            // Gas bond closed: combined fuel+clad axial equilibrium.
            for (int i = 0; i < nc; ++i)
                sum += m_geom.area0[nf + i] * s.sigz[i];
            sum += pcool * PI * m_geom.rco0 * m_geom.rco0;
        }
        r[n++] = sum;
    }

    // --- Fuel BC-2: inner radial boundary ---
    if (m_geom.rfi0 == 0.0) {
        // Solid pellet (no central hole): symmetry → sigfr[0] = sigfh[0]
        r[n++] = s.sigfr[0] - s.sigfh[0];
    } else {
        // Central hole: gas/bond pressure acts → sigfr[0] = -gpres
        r[n++] = s.sigfr[0] + gpres;
    }

    // --- Fuel BC-3: outer radial / gap interface ---
    if (s.gapOpen) {
        // Gas bond open: outer fuel surface exposed to fill gas → sigfr[nf-1] = -gpres
        r[n++] = s.sigfr[nf-1] + gpres;
    } else {
        // Gas bond closed: preserve axial strain offset frozen at gap closure.
        // axialOffset = efz - ez recorded at the closure event (set by the solver).
        // Algebraically equivalent to the legacy defz=dez rate constraint.
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
    // ec[i] is the hoop creep strain (ODE state, 0 for FRED-ROD).
    // Volume-preserving decomposition: ec_hoop = +ec, ec_radial = ec_axial = -ec/2.
    // Hoop:   E * (eh - et - ec) = sigh - nu*(sigr + sigz)
    for (int i = 0; i < nc; ++i) {
        double eel = s.eh[i] - s.et[i] - s.ec[i];
        r[n++] = Ec[i]*eel - s.sigh[i] + nuC[i]*(s.sigr[i] + s.sigz[i]);
    }
    // Radial: E * (er - et + ec/2) = sigr - nu*(sigh + sigz)
    for (int i = 0; i < nc; ++i) {
        double eel = s.er[i] - s.et[i] + 0.5*s.ec[i];
        r[n++] = Ec[i]*eel - s.sigr[i] + nuC[i]*(s.sigh[i] + s.sigz[i]);
    }
    // Axial:  E * (ez - et + ec/2) = sigz - nu*(sigh + sigr)
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
    if (s.gapOpen) {
        // Gap open: clad axial equilibrium (pressure difference between coolant
        //           and fill gas / bond pressure acts on the cladding annular ends)
        double sum = 0.0;
        for (int i = 0; i < nc; ++i)
            sum += m_geom.area0[nf + i] * s.sigz[i];
        sum += pcool * PI * m_geom.rco0 * m_geom.rco0
               - gpres * PI * m_geom.rci0 * m_geom.rci0;
        r[n++] = sum;
    } else {
        // Gap closed: radial stress continuity at fuel–clad interface.
        r[n++] = s.sigfr[nf-1] - s.sigr[0];
    }

    // --- Clad BC-2: inner radial surface ---
    if (s.gapOpen) {
        // Gap open: bond/fill-gas pressure on inner clad surface
        r[n++] = s.sigr[0] + gpres;
    } else {
        // Gap closed: geometric contact constraint (gap = roughness)
        r[n++] = s.gap - (m_geom.ruff + m_geom.rufc);
    }

    // --- Clad BC-3: outer radial surface (coolant pressure) ---
    r[n++] = s.sigr[nc-1] + pcool;

    assert(n == neq());
}

} // namespace fred
