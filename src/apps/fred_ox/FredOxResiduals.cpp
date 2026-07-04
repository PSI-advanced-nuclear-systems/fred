#include "FredOxResiduals.hpp"
#include "FredOxIrradiationPhysics.hpp"
#include "fuelpelletmaterial/FredOxMOX.hpp"
#include "platform/CladdingMaterial.hpp"
#include "platform/Constants.hpp"
#include "platform/GapModel.hpp"
#include "platform/GapPressureModel.hpp"
#include <cassert>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <omp.h>

namespace {

double computeOxGapConductance(
    const fred::FredOxLayerState& s,
    double gpres,
    const fred::FuelRodGeometry& geom,
    const fred::FredOxGapMaterial& gap_ox,
    const fred::FuelPelletMaterial& fuel,
    const fred::CladdingMaterial& clad)
{
    const int nf = geom.nf;
    const double T_gap_avg = 0.5 * (s.T[nf-1] + s.T[nf]);
    const double k_gap = gap_ox.gapConductivity(T_gap_avg);

    double h_gap = fred::computeGapMediumConductance(
        s.T[nf-1], s.T[nf],
        s.gap, gpres,
        k_gap,
        geom.ruff, geom.rufc,
        gap_ox.relocationFraction());
    h_gap += fred::computeRadiationConductance(
        s.T[nf-1], s.T[nf],
        geom.rfo0, geom.rci0);

    if (!s.gapOpen) {
        h_gap += fred::computeContactConductance(
            s.T[nf-1], s.T[nf],
            0.0, s.pfc, gpres,
            geom.ruff, geom.rufc,
            fuel, clad);
    }
    return h_gap;
}

}

namespace fred {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
FredOxResiduals::FredOxResiduals(const FuelRodGeometry& geom,
                                  FredOxMOX&             mox,
                                  const CladdingMaterial& clad,
                                  FredOxGapMaterial&      gap_ox,
                                  double coolant_pressure_MPa,
                                  double rof0_kg_m3,
                                  double pu_content,
                                  double sto0,
                                  double swelling_multiplier)
    : RodResiduals(geom, mox, clad, gap_ox, coolant_pressure_MPa,
                   /*neq_irr=*/1 + geom.nf + geom.nc),
      m_mox     (mox),
      m_gap_ox  (gap_ox),
      m_rof0    (rof0_kg_m3),
      m_pu      (pu_content),
      m_sto0    (sto0),
      m_fswelmlt(swelling_multiplier),
      m_state   (geom.nf, geom.nc, geom.nz)
{
    // As-fabricated free volume (gap annulus + central hole + plenum)
    double h_total = 0.0;
    for (int j = 0; j < geom.nz; ++j) h_total += geom.dz0[j];
    m_vfree0 = PI * (geom.rci0*geom.rci0 - geom.rfo0*geom.rfo0 + geom.rfi0*geom.rfi0)
               * h_total + geom.vgp;
}

// ---------------------------------------------------------------------------
// Gap state management
// ---------------------------------------------------------------------------
void FredOxResiduals::setGapOpen(int layer, bool open) {
    m_gapMgr.setGapOpen(layer, open);
    m_gapMgr.applyToLayer(layer, m_state.layers[layer]);
}

void FredOxResiduals::onGapClosed(int layer) {
    m_gapMgr.applyGapClosed(layer, m_state.layers[layer].efz, m_state.layers[layer].ez);
    m_gapMgr.applyToLayer(layer, m_state.layers[layer]);
}

// ---------------------------------------------------------------------------
// Gas pressure (ideal gas + fission gas)
// ---------------------------------------------------------------------------
double FredOxResiduals::computeGasPressure(const FredOxRodState& state) const {
    const double vt = computeGasBondVT(m_geom, state.layers, m_Tplenum);
    if (vt <= 0.0) return m_gpres0;
    return (m_mu0 + state.fgrel) * R_GAS / vt * 1.0e-6;
}

// ---------------------------------------------------------------------------
// fillIdArray
// ---------------------------------------------------------------------------
void FredOxResiduals::fillIdArray(double* id) const {
    const int nf = m_geom.nf, nc = m_geom.nc, nz = m_geom.nz;
    // Global: fggen(1), fgrel(1), gpres(0)
    id[0] = 1.0; id[1] = 1.0; id[2] = 0.0;
    for (int j = 0; j < nz; ++j) {
        double* idj = id + 3 + j * m_neq_j;
        // Thermal AEs
        idj[0] = 0.0; idj[1] = 0.0; idj[2] = 0.0; idj[3] = 0.0;
        // Fuel T ODEs
        for (int i = 0; i < nf; ++i)    idj[4+i] = 1.0;
        // Inner clad T ODEs
        for (int i = 0; i < nc-1; ++i)  idj[4+nf+i] = 1.0;
        // Outer clad T AE
        idj[4+nf+nc-1] = 0.0;
        // Burnup ODE
        idj[m_neq_th] = 1.0;
        // Fuel swelling ODEs
        for (int i = 0; i < nf; ++i) idj[offEfs(i)] = 1.0;
        // Cladding creep ODEs
        for (int i = 0; i < nc; ++i) idj[offEc (i)] = 1.0;
        // Mechanics: all AEs
        for (int k = offMech(); k < m_neq_j; ++k) idj[k] = 0.0;
    }
}

// ---------------------------------------------------------------------------
// pack / unpack
// ---------------------------------------------------------------------------
void FredOxResiduals::pack(const FredOxRodState& state, double* y) const {
    const int nf = m_geom.nf, nc = m_geom.nc;
    // Global block
    y[0] = state.fggen;
    y[1] = state.fgrel;
    y[2] = state.gpres;

    for (int j = 0; j < m_geom.nz; ++j) {
        double* yj = y + 3 + j * m_neq_j;
        const auto& s = state.layers[j];

        // Thermal block
        yj[0] = s.qqv;
        yj[1] = s.gap;
        yj[2] = s.pfc;
        yj[3] = s.hgap;
        for (int i = 0; i < nf+nc; ++i) yj[4+i] = s.T[i];

        // Irradiation block
        yj[m_neq_th] = s.bup * BUP_SCALE;
        for (int i = 0; i < nf; ++i) yj[offEfs(i)] = s.efs[i];
        for (int i = 0; i < nc; ++i) yj[offEc (i)] = s.ec[i];

        // Mechanics block (same layout as FredRodResiduals)
        for (int i = 0; i < nf; ++i) yj[offEft(i)]   = s.eft[i];
        for (int i = 0; i < nf; ++i) yj[offEfh(i)]   = s.efh[i];
        for (int i = 0; i < nf; ++i) yj[offEfr(i)]   = s.efr[i];
        yj[offEfz()]                                 = s.efz;
        for (int i = 0; i < nf; ++i) yj[offSigfh(i)] = s.sigfh[i];
        for (int i = 0; i < nf; ++i) yj[offSigfr(i)] = s.sigfr[i];
        for (int i = 0; i < nf; ++i) yj[offSigfz(i)] = s.sigfz[i];
        for (int i = 0; i < nc; ++i) yj[offEt(i)]    = s.et[i];
        for (int i = 0; i < nc; ++i) yj[offEh(i)]    = s.eh[i];
        for (int i = 0; i < nc; ++i) yj[offEr(i)]    = s.er[i];
        yj[offEz()]                                  = s.ez;
        for (int i = 0; i < nc; ++i) yj[offSigh(i)]  = s.sigh[i];
        for (int i = 0; i < nc; ++i) yj[offSigr(i)]  = s.sigr[i];
        for (int i = 0; i < nc; ++i) yj[offSigz(i)]  = s.sigz[i];
    }
}

void FredOxResiduals::unpack(const double* y, const double* yp, FredOxRodState& state) const {
    const int nf = m_geom.nf, nc = m_geom.nc;
    state.fggen  = y[0];
    state.fgrel  = std::max(0.0, y[1]);
    state.gpres  = y[2];
    state.dfggen = yp ? yp[0] : 0.0;
    state.dfgrel = yp ? yp[1] : 0.0;

    for (int j = 0; j < m_geom.nz; ++j) {
        const double* yj  = y  + 3 + j * m_neq_j;
        const double* ypj = yp ? yp + 3 + j * m_neq_j : nullptr;
        auto& s = state.layers[j];

        s.qqv  = yj[0];
        s.gap  = yj[1];
        s.pfc  = yj[2];
        s.hgap = yj[3];
        for (int i = 0; i < nf+nc; ++i) s.T[i] = yj[4+i];
        if (ypj) for (int i = 0; i < nf+nc; ++i) s.dT[i] = ypj[4+i];

        // Irradiation
        s.bup  = yj[m_neq_th] / BUP_SCALE;
        s.dbup = ypj ? ypj[m_neq_th] / BUP_SCALE : 0.0;
        for (int i = 0; i < nf; ++i) {
            s.efs[i]  = yj[offEfs(i)];
            s.defs[i] = ypj ? ypj[offEfs(i)] : 0.0;
        }
        for (int i = 0; i < nc; ++i) {
            s.ec[i]   = yj[offEc(i)];
            s.dec[i]  = ypj ? ypj[offEc(i)] : 0.0;
        }

        // Mechanics
        for (int i = 0; i < nf; ++i) s.eft[i]   = yj[offEft(i)];
        for (int i = 0; i < nf; ++i) s.efh[i]   = yj[offEfh(i)];
        for (int i = 0; i < nf; ++i) s.efr[i]   = yj[offEfr(i)];
        s.efz                                   = yj[offEfz()];
        for (int i = 0; i < nf; ++i) s.sigfh[i] = yj[offSigfh(i)];
        for (int i = 0; i < nf; ++i) s.sigfr[i] = yj[offSigfr(i)];
        for (int i = 0; i < nf; ++i) s.sigfz[i] = yj[offSigfz(i)];
        for (int i = 0; i < nc; ++i) s.et[i]    = yj[offEt(i)];
        for (int i = 0; i < nc; ++i) s.eh[i]    = yj[offEh(i)];
        for (int i = 0; i < nc; ++i) s.er[i]    = yj[offEr(i)];
        s.ez                                    = yj[offEz()];
        for (int i = 0; i < nc; ++i) s.sigh[i]  = yj[offSigh(i)];
        for (int i = 0; i < nc; ++i) s.sigr[i]  = yj[offSigr(i)];
        for (int i = 0; i < nc; ++i) s.sigz[i]  = yj[offSigz(i)];

        m_gapMgr.applyToLayer(j, s);

        // Update deformed geometry (needed for gas pressure and gap conductance)
        double dzf = 0.0, dzc = 0.0;
        s.updateGeometry(m_geom.rad0, m_geom.drf0, m_geom.drc0, m_geom.dz0[j], dzf, dzc);

        // Effective stress for fuel and clad
        for (int i = 0; i < nf; ++i) {
            double h = s.sigfh[i], r = s.sigfr[i], z = s.sigfz[i];
            s.sigf[i] = std::sqrt((h-r)*(h-r) + (h-z)*(h-z) + (z-r)*(z-r)) / std::sqrt(2.0);
        }
        for (int i = 0; i < nc; ++i) {
            double h = s.sigh[i], r = s.sigr[i], z = s.sigz[i];
            s.sig[i] = std::sqrt((h-r)*(h-r) + (h-z)*(h-z) + (z-r)*(z-r)) / std::sqrt(2.0);
        }
    }
}

// ---------------------------------------------------------------------------
// Global residuals (fggen, fgrel, gpres)
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// RodResiduals virtual hook implementations
// ---------------------------------------------------------------------------
void FredOxResiduals::unpackState(const double* y, const double* yp) const {
    unpack(y, yp, m_state);
    // Update gas inventory and plenum temperature on every call.
    m_gap_ox.setGasInventory(m_mu0);
    if (m_plenumTFn) m_Tplenum = m_plenumTFn(0.0); // updated per-t below in computeGlobalResiduals
}

void FredOxResiduals::prepareLayer(int j) const {
    m_state.layers[j].gapOpen = m_gapMgr.isGapOpen(j);
}

void FredOxResiduals::prepareThreadLocalState(int nthreads) const {
    const int n = nthreads > 0 ? nthreads : 1;
    if ((int)m_mox_tls.size() == n) return;
    m_mox_tls.assign(n, m_mox);
    m_gap_ox_tls.assign(n, m_gap_ox);
    m_heat_tls.clear();
    m_heat_tls.reserve(n);
    for (int t = 0; t < n; ++t)
        m_heat_tls.push_back(std::make_unique<HeatConduction>(m_geom, m_mox_tls[t], m_clad, m_gap_ox_tls[t]));
}

void FredOxResiduals::computeGlobalResiduals(double t, const double* yp, double* r) const
{
    // Refresh plenum temperature for this time step.
    if (m_plenumTFn) m_Tplenum = m_plenumTFn(t);

    const int nf = m_geom.nf;
    double fgGrate = 0.0, fgRrate = 0.0;

    // Hot start (isIrradiationOn()==false): freeze fggen/fgrel by leaving
    // both rates at 0, forcing yp[0]/yp[1]=0 below.
    if (isIrradiationOn()) {
        for (int j = 0; j < m_geom.nz; ++j) {
            const auto& s = m_state.layers[j];
            const double vol       = m_geom.azf0 * m_geom.dz0[j];
            const double ql_Wm     = s.qqv * m_geom.azf0;
            const double rate_fiss = fissionRate(s.qqv) * vol;
            const double gen_rate  = 0.25 * rate_fiss / AVOGADRO_OX;

            fgGrate += gen_rate;
            fgRrate += gen_rate * fgReleaseRate(s.bup, s.T.data(), m_geom.area0.data(),
                                                 nf, m_pu, ql_Wm);
        }
        if (fgGrate < 0.0) { fgGrate = 0.0; fgRrate = 0.0; }
    }

    r[0] = yp[0] - fgGrate;
    r[1] = yp[1] - fgRrate;
    r[2] = m_state.gpres - computeGasPressure(m_state);
}

// ---------------------------------------------------------------------------
// Thermal residuals for layer j
// ---------------------------------------------------------------------------
void FredOxResiduals::computeThermalResiduals(int j, double t,
                                               const double* /*yj*/,
                                               const double* /*ypj*/,
                                               double* rj) const
{
    const auto& s    = m_state.layers[j];
    const double qqv    = layerPower(j, t);
    const double T_cool = layerCoolant(j, t);

    // Per-thread clones (see prepareThreadLocalState): FredOxMOX/
    // FredOxGapMaterial are mutated then immediately read back below, so
    // each OpenMP worker thread must use its own instance, not the shared
    // m_mox/m_gap_ox (which would otherwise be a cross-layer data race).
    const int tid = (m_mox_tls.size() > 1) ? omp_get_thread_num() : 0;
    FredOxMOX&         mox    = m_mox_tls[tid];
    FredOxGapMaterial& gap_ox = m_gap_ox_tls[tid];

    mox.setBurnup(s.bup);
    mox.setStoichiometry(m_sto0);

    gap_ox.setFissionGasRelease(m_state.fgrel);
    gap_ox.setGasPressure(m_state.gpres);
    gap_ox.setBurnup(s.bup);
    gap_ox.setLinearPower(s.qqv * PI * m_geom.azf0);

    std::vector<double>& r_th = threadLocalThermalBuffer();
    m_heat_tls[tid]->computeResiduals(s, T_cool, m_state.gpres, qqv, r_th);
    r_th[3] = s.hgap - computeOxGapConductance(s, m_state.gpres, m_geom, gap_ox,
                                                 mox, m_clad);
    std::copy(r_th.begin(), r_th.end(), rj);
}

// ---------------------------------------------------------------------------
// Irradiation residuals (burnup + swelling) for layer j
// ---------------------------------------------------------------------------
void FredOxResiduals::computeIrradiationResiduals(int j,
                                                    const double* ypj,
                                                    double* rj) const
{
    const auto& s = m_state.layers[j];

    // Burnup ODE: d(bup_y)/dt = qqv/rof0
    rj[m_neq_th] = ypj[m_neq_th] - s.qqv / m_rof0;

    // Fuel swelling ODEs: d(efs[i])/dt = swelling_rate * multiplier
    for (int i = 0; i < m_geom.nf; ++i) {
        const double swel_rate = fuelSwellingRate(s.T[i], s.bup, s.qqv, m_rof0) * m_fswelmlt;
        rj[offEfs(i)] = ypj[offEfs(i)] - swel_rate;
    }

    // Cladding creep ODEs: d(ec[i])/dt = creepRate(T[nf+i], sigh[i])
    // ec[i] is the hoop creep strain; radial/axial creep follow volume-preserving
    // decomposition (handled in StressStrain::computeResiduals).
    for (int i = 0; i < m_geom.nc; ++i) {
        const double rate = m_clad.creepRate(s.T[m_geom.nf + i], s.sigh[i]);
        rj[offEc(i)] = ypj[offEc(i)] - rate;
    }
}

// ---------------------------------------------------------------------------
// initAlgebraicState
// ---------------------------------------------------------------------------
void FredOxResiduals::initAlgebraicState(FredOxRodState& state) const {
    const int nf = m_geom.nf, nc = m_geom.nc;
    state.gpres  = m_gpres0;
    state.mu0    = m_mu0;
    state.fggen  = 0.0;
    state.fgrel  = 0.0;

    for (int j = 0; j < m_geom.nz; ++j) {
        auto& s = state.layers[j];
        s.qqv    = (!m_layerPowerFns.empty() && j < (int)m_layerPowerFns.size())
                   ? m_layerPowerFns[j](0.0)
                   : (m_powerFn ? m_powerFn(0.0) : 0.0);
        s.bup    = 0.0;
        s.gap    = m_geom.rci0 - m_geom.rfo0;
        s.pfc    = 0.0;
        m_gapMgr.applyToLayer(j, s);
        for (int i = 0; i < nf; ++i) s.efs[i] = 0.0;

        double Tc0 = (!m_layerCoolantFns.empty() && j < (int)m_layerCoolantFns.size())
                     ? m_layerCoolantFns[j](0.0)
                     : (m_coolantTFn ? m_coolantTFn(0.0) : m_T_init);
        s.T[nf+nc-1] = Tc0;

        // Lamé initial elastic solution (same as FredRodResiduals)
        m_mox.setBurnup(0.0);
        m_gap_ox.setFissionGasRelease(0.0);
        m_gap_ox.setGasPressure(m_gpres0);
        m_gap_ox.setLinearPower(0.0);

        s.hgap = computeOxGapConductance(s, m_gpres0, m_geom, m_gap_ox, m_mox, m_clad);

        if (m_mech_on && m_gapMgr.isGapOpen(j)) {
            for (int i = 0; i < nf; ++i) {
                double Ef  = m_mox.youngsModulus(s.T[i], m_mox.referenceDensity());
                double nuF = m_mox.poissonRatio();
                s.eft[i]   = m_mox.thermalExpansionStrain(s.T[i]);
                s.sigfr[i] = -m_gpres0;
                s.sigfh[i] = -m_gpres0;
                s.sigfz[i] = -m_gpres0;
                double e_mech = -m_gpres0 * (1.0 - 2.0*nuF) / Ef;
                s.efh[i] = s.eft[i] + e_mech;
                s.efr[i] = s.eft[i] + e_mech;
                s.efz    = s.eft[i] + e_mech;
            }
            double a = m_geom.rci0, b = m_geom.rco0;
            double a2 = a*a, b2 = b*b, denom = b2 - a2;
            double pcool_0 = m_pcoolFn ? m_pcoolFn(0.0) : m_pcool;
            double A_c = (m_gpres0*a2 - pcool_0*b2) / denom;
            double C_c = (pcool_0 - m_gpres0)*a2*b2 / denom;
            for (int i = 0; i < nc; ++i) {
                double ri  = m_geom.rad0[nf+i], ri2 = ri*ri;
                double Ec  = m_clad.youngsModulus(s.T[nf+i]);
                double nuC = m_clad.poissonRatio();
                s.et[i]    = m_clad.thermalExpansionStrain(s.T[nf+i]);
                s.sigr[i]  = A_c + C_c/ri2;
                s.sigh[i]  = A_c - C_c/ri2;
                s.sigz[i]  = A_c;
                s.eh[i]    = s.et[i] + (s.sigh[i] - nuC*(s.sigr[i]+s.sigz[i])) / Ec;
                s.er[i]    = s.et[i] + (s.sigr[i] - nuC*(s.sigh[i]+s.sigz[i])) / Ec;
                s.ez       = s.et[i] + (s.sigz[i] - nuC*(s.sigh[i]+s.sigr[i])) / Ec;
            }
        }

        // Update deformed geometry
        double dzf = 0.0, dzc = 0.0;
        s.updateGeometry(m_geom.rad0, m_geom.drf0, m_geom.drc0, m_geom.dz0[j], dzf, dzc);
    }
}

void FredOxResiduals::solveMechanicalIC(FredOxRodState& state, double t,
                                         double rtol, double atol) const {
    // Delegate to a simple Newton solve on the mechanical equations.
    // For now, the Lamé seed from initAlgebraicState is a good enough starting point.
    (void)state; (void)t; (void)rtol; (void)atol;
}

} // namespace fred
