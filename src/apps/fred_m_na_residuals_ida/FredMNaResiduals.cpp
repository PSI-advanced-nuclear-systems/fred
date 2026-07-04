#include "FredMNaResiduals.hpp"
#include "FredMNaIrradiationPhysics.hpp"
#include "platform/CladdingMaterial.hpp"
#include "platform/Constants.hpp"
#include "platform/GapModel.hpp"
#include "platform/GapMaterial.hpp"
#include "platform/GapPressureModel.hpp"
#include <nvector/nvector_serial.h>
#include <cassert>
#include <cmath>
#include <algorithm>

namespace {

double computeMNaGapConductance(
    const fred::FredMNaLayerState& s,
    double gpres,
    const fred::FuelRodGeometry& geom,
    const fred::GapMaterial& gap_mat,
    const fred::FuelPelletMaterial& fuel,
    const fred::CladdingMaterial& clad)
{
    const int nf = geom.nf;
    const double T_gap_avg = 0.5 * (s.T[nf-1] + s.T[nf]);
    const double k_gap = gap_mat.gapConductivity(T_gap_avg);
    const double rough = std::sqrt(geom.ruff * geom.ruff + geom.rufc * geom.rufc);
    const double gap_eff = std::max(s.gap, rough);
    double h_gap = k_gap / gap_eff;

    if (!s.gapOpen) {
        h_gap += fred::computeContactConductance(
            s.T[nf-1], s.T[nf],
            0.0, s.pfc, gpres,
            geom.ruff, geom.rufc,
            fuel, clad);
    }
    h_gap = gap_mat.clampConductance(h_gap, !s.gapOpen);
    return h_gap;
}

} // namespace

namespace fred {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
FredMNaResiduals::FredMNaResiduals(const FuelRodGeometry& geom,
                                    UPuZr&                 fuel,
                                    const HT9&             clad,
                                    GapMaterial&           gap_mat,
                                    double                 coolant_pressure_MPa)
    : RodResiduals(geom, fuel, clad, gap_mat, coolant_pressure_MPa,
                   /*neq_irr=*/1 + geom.nf + geom.nc),
      m_fuel   (fuel),
      m_ht9    (clad),
      m_gap_mat(gap_mat),
      m_state  (geom.nf, geom.nc, geom.nz),
      m_mech_mna(geom, fuel, clad),
      m_r_mech_mna(m_mech_mna.neq())
{
    for (int j = 0; j < geom.nz; ++j) m_h_total += geom.dz0[j];
    // Liquid bond: fill gas confined to plenum only.
    m_vfree0 = geom.vgp;
}

// ---------------------------------------------------------------------------
// Gap state management
// ---------------------------------------------------------------------------
void FredMNaResiduals::setGapOpen(int layer, bool open) {
    m_gapMgr.setGapOpen(layer, open);
    m_gapMgr.applyToLayer(layer, m_state.layers[layer]);
}

void FredMNaResiduals::onGapClosed(int layer) {
    // Use setGapOpen(false) rather than applyGapClosed to keep axialOffset = 0.
    // FredMNa uses efz = ez at closure, so StressStrain BC becomes efz - ez - 0 = 0.
    m_gapMgr.setGapOpen(layer, false);
    m_gapMgr.applyToLayer(layer, m_state.layers[layer]);
}

// ---------------------------------------------------------------------------
// Gas pressure — Precal.for's "upuzr" (sodium/liquid-bond metal fuel) branch.
//
// Unlike the general oxide ideal-gas-law model (fresh mu/V/T each step), the
// legacy metal-fuel formula is additive on top of the fixed backfill pressure
// gpres0: as fuel swells and the gap shrinks, sodium bond fluid is displaced
// out of the annular gap and up into the plenum, shrinking the plenum's own
// free (gas) volume (delfv, using the fuel radius capped at the clad ID once
// in contact, and the clad OD's own growth). Gas trapped in GRSIS open-bubble
// porosity within the fuel (swopen = sw_v3) adds its own V/T term at the
// local per-layer average fuel temperature (vt_axn). No mu0/fggen term here —
// this matches Precal.for exactly, where gpres = gpres0 + fgrel*R/(vpl/tple +
// vt_axn).
// ---------------------------------------------------------------------------
double FredMNaResiduals::computeGasPressure(const FredMNaRodState& state) const {
    const int nf = m_geom.nf;
    double delfv = 0.0, vt_axn = 0.0;
    for (int j = 0; j < m_geom.nz; ++j) {
        const auto& s    = state.layers[j];
        const double dzf = m_geom.dz0[j] * (1.0 + s.efz);
        const double dzc = m_geom.dz0[j] * (1.0 + s.ez);
        const double radius_fo = std::min(s.rfo, s.rci);

        delfv += PI*radius_fo*radius_fo*dzf - PI*m_geom.rfo0*m_geom.rfo0*m_geom.dz0[j]
               - PI*s.rco*s.rco*dzc          + PI*m_geom.rco0*m_geom.rco0*m_geom.dz0[j];

        double T_fave = 0.0;
        for (int i = 0; i < nf; ++i) T_fave += s.T[i];
        T_fave = std::max(T_fave / nf, 1.0);

        double v_void = 0.0;
        for (int i = 0; i < nf-1; ++i)
            v_void += s.nodes[i].grsis.swopen
                    * PI * (s.rad[i+1]*s.rad[i+1] - s.rad[i]*s.rad[i]) * dzf;
        vt_axn += v_void / T_fave;
    }
    const double vpl    = m_geom.vgp - delfv;
    const double vt_eff = vpl / std::max(m_Tplenum, 1.0) + vt_axn;
    if (vt_eff <= 0.0) return m_gpres0;
    return m_gpres0 + state.fgrel * R_GAS / vt_eff * 1.0e-6;
}

// ---------------------------------------------------------------------------
// syncAuxLayerState — see header comment.
// ---------------------------------------------------------------------------
void FredMNaResiduals::syncAuxLayerState(int j, const FredMNaLayerState& src) {
    auto& dst = m_state.layers[j];
    dst.k_irr_factor  = src.k_irr_factor;
    dst.flag          = src.flag;
    dst.bup_FIMA      = src.bup_FIMA;
    dst.buhard_FIMA   = src.buhard_FIMA;
    dst.xwast         = src.xwast;
    dst.clfuel        = src.clfuel;
    dst.ntot          = src.ntot;
    dst.defs_grsis    = src.defs_grsis;
    dst.efsz          = src.efsz;
    dst.efsh          = src.efsh;
    dst.efsr          = src.efsr;
    for (int i = 0; i < m_geom.nf; ++i) {
        auto&       nd  = dst.nodes[i];
        const auto& snd = src.nodes[i];
        nd.zr_wf   = snd.zr_wf;   nd.pu_wf   = snd.pu_wf;   nd.ur_wf = snd.ur_wf;
        nd.zr_at   = snd.zr_at;   nd.pu_at   = snd.pu_at;   nd.ur_at = snd.ur_at;
        nd.c_zr    = snd.c_zr;
        nd.phase   = snd.phase;   nd.pfrac   = snd.pfrac;   nd.psod  = snd.psod;
        nd.poros_tot = snd.poros_tot;
        nd.poros_gas = snd.poros_gas;
        nd.grsis     = snd.grsis;
    }
}

// ---------------------------------------------------------------------------
// fillIdArray
// ---------------------------------------------------------------------------
void FredMNaResiduals::fillIdArray(double* id) const {
    const int nf = m_geom.nf, nc = m_geom.nc, nz = m_geom.nz;
    id[0] = 1.0; id[1] = 0.0; id[2] = 0.0; // fggen ODE, fgrel AE (GRSIS), gpres AE
    for (int j = 0; j < nz; ++j) {
        double* idj = id + 3 + j * m_neq_j;
        idj[0] = 0.0; idj[1] = 0.0; idj[2] = 0.0; idj[3] = 0.0; // thermal AEs
        for (int i = 0; i < nf; ++i)   idj[4+i]    = 1.0; // fuel T ODEs
        for (int i = 0; i < nc-1; ++i) idj[4+nf+i] = 1.0; // inner clad T ODEs
        idj[4+nf+nc-1] = 0.0;                              // outer clad T AE
        idj[m_neq_th] = 1.0;                               // bup_y ODE
        for (int i = 0; i < nf; ++i) idj[offEfs(i)] = 1.0; // efs ODEs
        for (int i = 0; i < nc; ++i) idj[offEc (i)] = 1.0; // ec ODEs
        for (int k = offMech(); k < m_neq_j; ++k) idj[k] = 0.0; // mechanics AEs
    }
}

// ---------------------------------------------------------------------------
// fillFdScale — see header doc comment. One layer's worth (size m_neq_j);
// caller (FredMNaSolver::newtonSolveLayer) indexes into this per-layer, and
// appends its own scale entry for the jointly-solved gpres unknown.
// ---------------------------------------------------------------------------
void FredMNaResiduals::fillFdScale(double* scale) const {
    const int nf = m_geom.nf, nc = m_geom.nc;
    scale[0] = 1.0e9;   // qqv [W/m3] — peak power density order of magnitude
    scale[1] = 1.0e-4;  // gap [m] — typical as-fabricated gap, ~100 um
    scale[2] = 10.0;    // pfc [MPa] — contact pressure order of magnitude
    scale[3] = 1.0e4;   // hgap [W/m2K] — gap conductance order of magnitude
    for (int i = 0; i < nf+nc; ++i) scale[4+i] = 1.0e3;   // T [K]
    scale[m_neq_th] = 1.0e9;                              // bup_y = bup*MNA_BUP_SCALE
    for (int i = 0; i < nf; ++i) scale[offEfs(i)] = 1.0e-2; // efs [strain]
    for (int i = 0; i < nc; ++i) scale[offEc (i)] = 1.0e-3; // ec  [creep strain]

    for (int i = 0; i < nf; ++i) scale[offEft(i)]   = 1.0e-2; // fuel strains
    for (int i = 0; i < nf; ++i) scale[offEfh(i)]   = 1.0e-2;
    for (int i = 0; i < nf; ++i) scale[offEfr(i)]   = 1.0e-2;
    scale[offEfz()]                                   = 1.0e-2;
    for (int i = 0; i < nf; ++i) scale[offSigfh(i)] = 10.0;   // fuel stresses [MPa]
    for (int i = 0; i < nf; ++i) scale[offSigfr(i)] = 10.0;
    for (int i = 0; i < nf; ++i) scale[offSigfz(i)] = 10.0;
    for (int i = 0; i < nc; ++i) scale[offEt(i)]    = 1.0e-2; // clad strains
    for (int i = 0; i < nc; ++i) scale[offEh(i)]    = 1.0e-2;
    for (int i = 0; i < nc; ++i) scale[offEr(i)]    = 1.0e-2;
    scale[offEz()]                                    = 1.0e-2;
    for (int i = 0; i < nc; ++i) scale[offSigh(i)]  = 10.0;   // clad stresses [MPa]
    for (int i = 0; i < nc; ++i) scale[offSigr(i)]  = 10.0;
    for (int i = 0; i < nc; ++i) scale[offSigz(i)]  = 10.0;
}

// ---------------------------------------------------------------------------
// pack / unpack
// ---------------------------------------------------------------------------
void FredMNaResiduals::pack(const FredMNaRodState& state, double* y) const {
    const int nf = m_geom.nf, nc = m_geom.nc;
    y[0] = state.fggen;
    y[1] = state.fgrel;
    y[2] = state.gpres;

    for (int j = 0; j < m_geom.nz; ++j) {
        double* yj = y + 3 + j * m_neq_j;
        const auto& s = state.layers[j];

        yj[0] = s.qqv; yj[1] = s.gap; yj[2] = s.pfc; yj[3] = s.hgap;
        for (int i = 0; i < nf+nc; ++i) yj[4+i] = s.T[i];

        yj[m_neq_th] = s.bup * MNA_BUP_SCALE;
        for (int i = 0; i < nf; ++i) yj[offEfs(i)] = s.efs[i];
        for (int i = 0; i < nc; ++i) yj[offEc (i)] = s.ec[i];

        for (int i = 0; i < nf; ++i) yj[offEft(i)]   = s.eft[i];
        for (int i = 0; i < nf; ++i) yj[offEfh(i)]   = s.efh[i];
        for (int i = 0; i < nf; ++i) yj[offEfr(i)]   = s.efr[i];
        yj[offEfz()]                                   = s.efz;
        for (int i = 0; i < nf; ++i) yj[offSigfh(i)] = s.sigfh[i];
        for (int i = 0; i < nf; ++i) yj[offSigfr(i)] = s.sigfr[i];
        for (int i = 0; i < nf; ++i) yj[offSigfz(i)] = s.sigfz[i];
        for (int i = 0; i < nc; ++i) yj[offEt(i)]    = s.et[i];
        for (int i = 0; i < nc; ++i) yj[offEh(i)]    = s.eh[i];
        for (int i = 0; i < nc; ++i) yj[offEr(i)]    = s.er[i];
        yj[offEz()]                                    = s.ez;
        for (int i = 0; i < nc; ++i) yj[offSigh(i)]  = s.sigh[i];
        for (int i = 0; i < nc; ++i) yj[offSigr(i)]  = s.sigr[i];
        for (int i = 0; i < nc; ++i) yj[offSigz(i)]  = s.sigz[i];
    }
}

void FredMNaResiduals::unpackLayer(int j, const double* yj, const double* ypj,
                                    FredMNaLayerState& s) const {
    const int nf = m_geom.nf, nc = m_geom.nc;

    s.qqv  = yj[0]; s.gap  = yj[1];
    s.pfc  = yj[2]; s.hgap = yj[3];
    for (int i = 0; i < nf+nc; ++i) s.T[i]  = yj[4+i];
    if (ypj) for (int i = 0; i < nf+nc; ++i) s.dT[i] = ypj[4+i];

    s.bup  = yj[m_neq_th] / MNA_BUP_SCALE;
    s.dbup = ypj ? ypj[m_neq_th] / MNA_BUP_SCALE : 0.0;
    for (int i = 0; i < nf; ++i) {
        s.efs[i]  = yj[offEfs(i)];
        s.defs[i] = ypj ? ypj[offEfs(i)] : 0.0;
    }
    for (int i = 0; i < nc; ++i) {
        s.ec[i]   = yj[offEc(i)];
        s.dec[i]  = ypj ? ypj[offEc(i)] : 0.0;
    }

    for (int i = 0; i < nf; ++i) s.eft[i]   = yj[offEft(i)];
    for (int i = 0; i < nf; ++i) s.efh[i]   = yj[offEfh(i)];
    for (int i = 0; i < nf; ++i) s.efr[i]   = yj[offEfr(i)];
    s.efz                                     = yj[offEfz()];
    for (int i = 0; i < nf; ++i) s.sigfh[i] = yj[offSigfh(i)];
    for (int i = 0; i < nf; ++i) s.sigfr[i] = yj[offSigfr(i)];
    for (int i = 0; i < nf; ++i) s.sigfz[i] = yj[offSigfz(i)];
    for (int i = 0; i < nc; ++i) s.et[i]    = yj[offEt(i)];
    for (int i = 0; i < nc; ++i) s.eh[i]    = yj[offEh(i)];
    for (int i = 0; i < nc; ++i) s.er[i]    = yj[offEr(i)];
    s.ez                                      = yj[offEz()];
    for (int i = 0; i < nc; ++i) s.sigh[i]  = yj[offSigh(i)];
    for (int i = 0; i < nc; ++i) s.sigr[i]  = yj[offSigr(i)];
    for (int i = 0; i < nc; ++i) s.sigz[i]  = yj[offSigz(i)];

    m_gapMgr.applyToLayer(j, s);

    double dzf = 0.0, dzc = 0.0;
    s.updateGeometry(m_geom.rad0, m_geom.drf0, m_geom.drc0, m_geom.dz0[j], dzf, dzc);

    for (int i = 0; i < nf; ++i) {
        double h = s.sigfh[i], r = s.sigfr[i], z = s.sigfz[i];
        s.sigf[i] = std::sqrt((h-r)*(h-r) + (h-z)*(h-z) + (z-r)*(z-r)) / std::sqrt(2.0);
    }
    for (int i = 0; i < nc; ++i) {
        double h = s.sigh[i], r = s.sigr[i], z = s.sigz[i];
        s.sig[i] = std::sqrt((h-r)*(h-r) + (h-z)*(h-z) + (z-r)*(z-r)) / std::sqrt(2.0);
    }
}

void FredMNaResiduals::unpack(const double* y, const double* yp,
                               FredMNaRodState& state) const {
    state.fggen  = y[0];
    // No clamp: fgrel is now an AE (r[1] = y[1] - m_fgrel_grsis) driven by GRSIS,
    // which is itself guaranteed non-negative. Clamping here would make the
    // residual constant (and the Jacobian row zero) whenever a Newton trial
    // iterate for y[1] goes negative, producing a singular iteration matrix.
    state.fgrel  = y[1];
    state.gpres  = y[2];
    state.dfggen = yp ? yp[0] : 0.0;
    state.dfgrel = yp ? yp[1] : 0.0;

    for (int j = 0; j < m_geom.nz; ++j) {
        const double* yj  = y  + 3 + j * m_neq_j;
        const double* ypj = yp ? yp + 3 + j * m_neq_j : nullptr;
        unpackLayer(j, yj, ypj, state.layers[j]);
    }
}

// ---------------------------------------------------------------------------
// initAlgebraicState
// ---------------------------------------------------------------------------
void FredMNaResiduals::initAlgebraicState(FredMNaRodState& state) const {
    const int nf = m_geom.nf, nc = m_geom.nc;
    state.fggen = 0.0;
    state.fgrel = 0.0;
    state.gpres = m_gpres0;

    for (int j = 0; j < m_geom.nz; ++j) {
        auto& s = state.layers[j];
        s.qqv  = m_layerPowerFns.empty() ? 0.0 : m_layerPowerFns[j](0.0);
        s.gap  = m_geom.rci0 - m_geom.rfo0;
        s.pfc  = 0.0;
        s.bup  = 0.0;
        s.bup_FIMA = 0.0;
        s.flag = s.gap > (m_geom.ruff + m_geom.rufc) ? "open" : "clos";
        for (int i = 0; i < nf; ++i) s.efs[i] = 0.0;

        for (int i = 0; i < nf+nc; ++i) s.T[i] = m_T_init;
        if (!m_layerCoolantFns.empty() && j < (int)m_layerCoolantFns.size())
            s.T[nf+nc-1] = m_layerCoolantFns[j](0.0);

        for (int i = 0; i < nf+nc; ++i) s.rad[i] = m_geom.rad0[i];
        s.rfi = m_geom.rfi0; s.rfo = m_geom.rfo0;
        s.rci = m_geom.rci0; s.rco = m_geom.rco0;

        const double pu = m_fuel.puContent();
        const double zr = m_fuel.zrContent();
        for (int i = 0; i < nf; ++i) {
            auto& nd = s.nodes[i];
            nd.pu_wf = pu; nd.zr_wf = zr; nd.ur_wf = 1.0 - pu - zr;
            nd.phase = "alpha"; nd.pfrac = 1.0; nd.psod = 0.0;
            nd.poros_tot = 0.0; nd.poros_gas = 0.0;
        }

        m_gapMgr.applyToLayer(j, s);
        s.hgap = computeMNaGapConductance(s, state.gpres, m_geom, m_gap_mat, m_fuel, m_clad);

        // Seed mechanical state at thermal-expansion equilibrium for IDACalcIC.
        {
            const double T0   = m_T_init;
            const double p0   = m_gpres0;
            const double eft0 = m_fuel.thermalExpansionStrain(T0);
            for (int i = 0; i < nf; ++i) {
                s.eft[i]   = eft0;
                s.efh[i]   = eft0;
                s.efr[i]   = eft0;
                s.sigfh[i] = -p0;
                s.sigfr[i] = -p0;
                s.sigfz[i] = -p0;
            }
            s.efz = eft0;

            const double et0 = m_clad.thermalExpansionStrain(T0);
            for (int i = 0; i < nc; ++i) {
                s.et[i]   = et0;
                s.eh[i]   = et0;
                s.er[i]   = et0;
                s.sigh[i] = -p0;
                s.sigr[i] = -p0;
                s.sigz[i] = -p0;
            }
            s.ez = et0;

            double dzf = 0.0, dzc = 0.0;
            s.updateGeometry(m_geom.rad0, m_geom.drf0, m_geom.drc0, m_geom.dz0[j], dzf, dzc);
        }
    }
}

// ---------------------------------------------------------------------------
// RodResiduals virtual hook: unpackState
// ---------------------------------------------------------------------------
void FredMNaResiduals::unpackState(const double* y, const double* yp) const {
    unpack(y, yp, m_state);
    if (m_plenumTFn) m_Tplenum = m_plenumTFn(0.0); // refined per-t in computeGlobalResiduals
}

// ---------------------------------------------------------------------------
// RodResiduals virtual hook: prepareLayer
// Per-node phase/infiltration/porosity update (bup_FIMA=0 during integration).
// ---------------------------------------------------------------------------
void FredMNaResiduals::prepareLayer(int j) const {
    auto& s = m_state.layers[j];
    const double pu = m_fuel.puContent();
    const double zr = m_fuel.zrContent();
    for (int i = 0; i < m_geom.nf; ++i) {
        auto& nd = s.nodes[i];
        auto ph = upuzrPhase(s.T[i], pu, zr);
        nd.phase = ph.phase; nd.pfrac = ph.pfrac;
        nd.psod = sodiumInfiltration(s.bup_FIMA, s.buhard_FIMA,
                                      m_geom.rad0[i], m_geom.rad0[m_geom.nf-1],
                                      nd.phase, s.flag);
        // Use GRSIS porosity if available (afterAcceptedStep has populated grsis.swtot),
        // otherwise fall back to the simple solid-FP swelling efs.
        if (nd.grsis.swtot > 0.0) {
            nd.poros_tot = nd.grsis.swtot;
            nd.poros_gas = nd.grsis.frtot;
        } else {
            nd.poros_tot = std::max(0.0, s.efs[i]);
            nd.poros_gas = 0.5 * nd.poros_tot;
        }
    }
    m_gapMgr.applyToLayer(j, s);
}

// ---------------------------------------------------------------------------
// RodResiduals virtual hook: computeGlobalResiduals
// ---------------------------------------------------------------------------
void FredMNaResiduals::computeGlobalResiduals(double t, const double* yp, double* r) const {
    // Metal-fuel plenum temperature = outlet (top-layer) coolant temperature,
    // recomputed every step (legacy Baseir.for: tple(l) = tcool(nzz(l),l) for
    // fmat=='upuzr'). An explicit override via setPlenumTemperature() still wins
    // if the caller supplied one.
    m_Tplenum = m_plenumTFn ? m_plenumTFn(t) : layerCoolant(m_geom.nz - 1, t);

    const double pu = m_fuel.puContent();
    const double zr = m_fuel.zrContent();

    const double molwe_HM = 1.0 / (pu / 0.244 + (1.0 - zr - pu) / 0.238029);
    const double bup_to_FIMA = (1.0e6 * 8.64e4 * molwe_HM)
                              / (3.204354e-11 * MNA_AVOGADRO);

    double bupave_FIMA = 0.0;
    for (int j = 0; j < m_geom.nz; ++j)
        bupave_FIMA += m_state.layers[j].bup * bup_to_FIMA;
    bupave_FIMA /= std::max(m_geom.nz, 1);

    double fgGrate = 0.0;
    for (int j = 0; j < m_geom.nz; ++j) {
        const auto& s = m_state.layers[j];
        const double vol      = m_geom.azf0 * m_geom.dz0[j];
        const double rate_fiss = s.qqv / (200.0 * 1.60214e-13);
        const double gen_rate  = 0.25 * rate_fiss * vol / MNA_AVOGADRO;
        fgGrate += gen_rate;
    }
    fgGrate = std::max(fgGrate, 0.0);

    r[0] = yp[0] - fgGrate;
    // fgrel is now an AE: driven by GRSIS cgrel set via setFgrelFromGrsis().
    r[1] = m_state.fgrel - m_fgrel_grsis;

    // gpres AE: ideal gas law, liquid-bond V/T (platform helper).
    r[2] = m_state.gpres - computeGasPressure(m_state);
}

// ---------------------------------------------------------------------------
// RodResiduals virtual hook: computeThermalResiduals
// ---------------------------------------------------------------------------
void FredMNaResiduals::computeThermalResiduals(int j, double t,
                                                const double* /*yj*/,
                                                const double* /*ypj*/,
                                                double* rj) const {
    const auto& s    = m_state.layers[j];
    const double qqv    = layerPower(j, t);
    const double T_cool = layerCoolant(j, t);

    const bool   prescribed = (m_subchannel == nullptr);
    const double h_cool_j   = prescribed ? 0.0 : layerHTC(j);
    m_heat.computeResiduals(s, T_cool, m_state.gpres, qqv, m_r_th, prescribed, h_cool_j);
    m_r_th[3] = s.hgap - computeMNaGapConductance(s, m_state.gpres, m_geom, m_gap_mat,
                                                    m_fuel, m_clad);
    std::copy(m_r_th.begin(), m_r_th.end(), rj);
}

// ---------------------------------------------------------------------------
// RodResiduals virtual hook: computeIrradiationResiduals
// ---------------------------------------------------------------------------
void FredMNaResiduals::computeIrradiationResiduals(int j,
                                                    const double* ypj,
                                                    double* rj) const {
    const auto& s   = m_state.layers[j];
    const double rof0 = m_fuel.referenceDensity();

    rj[m_neq_th] = ypj[m_neq_th] - s.qqv / rof0;

    const double pu    = m_fuel.puContent();
    const double zr    = m_fuel.zrContent();
    const double molwe = 1.0 / (pu / 0.244 + (1.0 - zr - pu) / 0.238029);
    const double hm    = rof0 * MNA_AVOGADRO / molwe;
    const double d_bup_FIMA_dt = s.qqv / (3.204354e-11 * hm);

    const double swel_rate_lin = 0.015 * 100.0 * d_bup_FIMA_dt / 3.0;

    // Once GRSIS has run for a node (nodes[i].grsis.swtot > 0), drive efs's rate
    // from GRSIS's total swelling (solid FP + gas bubbles) instead of the
    // solid-FP-only empirical formula, so gap closure (and thus cladding
    // wastage) is physically coupled to fission-gas swelling. Same fallback
    // pattern as prepareLayer's poros_tot (falls back to the empirical rate
    // when GRSIS is disabled or hasn't produced a value yet).
    for (int i = 0; i < m_geom.nf; ++i) {
        const double rate = (s.nodes[i].grsis.swtot > 0.0)
                           ? s.defs_grsis[i]
                           : swel_rate_lin;
        rj[offEfs(i)] = ypj[offEfs(i)] - rate;
    }

    // Cladding hoop creep ODEs: d(ec[i])/dt = creepRate(T_clad[i], sigh[i], qqv, t)
    // Clamp tensile stress to zero: pow(negative, non-integer) is NaN.
    for (int i = 0; i < m_geom.nc; ++i) {
        const double sigh_Pa = std::max(0.0, s.sigh[i]) * 1.0e6;
        const double rate    = m_ht9.creepRate(s.T[m_geom.nf + i],
                                                sigh_Pa,
                                                s.qqv,
                                                m_elapsed_time);
        rj[offEc(i)] = ypj[offEc(i)] - rate;
    }
}

// ---------------------------------------------------------------------------
// computeResidualsMNa / idaResidualMNa — FRED-M-Na's own IDA residual entry
// point (Decision D2/D8). Reimplements the same per-layer loop as
// RodResiduals::computeResiduals (platform/RodResiduals.cpp, untouched per
// Decision D1), but calls FredMNaStressStrain instead of the inherited
// (two-state) m_mech, so the three-state open/soft/clos gap-contact BCs and
// directional swelling partition are used for M-Na mechanics.
// ---------------------------------------------------------------------------
int FredMNaResiduals::computeResidualsMNa(double t, const double* y, const double* yp, double* r) const {
    unpackState(y, yp);

    const int ng = globalOffset(); // 3

    computeGlobalResiduals(t, yp, r);

    if (m_subchannel != nullptr) {
        std::vector<double> qqv_all(m_geom.nz);
        for (int j = 0; j < m_geom.nz; ++j)
            qqv_all[j] = layerPower(j, t);
        m_subchannel->updateCoolantField(t, qqv_all.data(), m_geom.dz0.data());
    }

    for (int j = 0; j < m_geom.nz; ++j) {
        const double* yj  = y  + ng + j * m_neq_j;
        const double* ypj = yp + ng + j * m_neq_j;
        double*       rj  = r  + ng + j * m_neq_j;
        computeLayerResiduals(j, t, yj, ypj, rj);
    }
    return 0;
}

int FredMNaResiduals::idaResidualMNa(double t, N_Vector y, N_Vector yp, N_Vector r, void* user_data) {
    auto* self = static_cast<FredMNaResiduals*>(user_data);
    return self->computeResidualsMNa(
        t,
        N_VGetArrayPointer(y),
        N_VGetArrayPointer(yp),
        N_VGetArrayPointer(r));
}

// ---------------------------------------------------------------------------
// computeLayerResiduals — one layer's thermal + irradiation + mechanics
// residual block only (O(1), not O(nz)). See header comment: used by the
// one-step backward-Euler integrator's per-layer local Newton solve, and by
// computeResidualsMNa above (single implementation, no duplication).
// Assumes m_state.gpres and all other layers are already current.
// ---------------------------------------------------------------------------
void FredMNaResiduals::computeLayerResiduals(int j, double t, const double* yj,
                                              const double* ypj, double* rj) const {
    unpackLayer(j, yj, ypj, m_state.layers[j]);
    prepareLayer(j);
    const double pcool = layerPcool(t);

    if (isHeatOn())
        computeThermalResiduals(j, t, yj, ypj, rj);
    else
        freezeThermalResiduals(j, yj, rj);

    computeIrradiationResiduals(j, ypj, rj);

    if (isMechOn()) {
        m_mech_mna.computeResiduals(m_state.layers[j], m_state.gpres, pcool, m_r_mech_mna);
        std::copy(m_r_mech_mna.begin(), m_r_mech_mna.end(), rj + m_neq_th + m_neq_irr);
    } else {
        // Pin each mechanical slot to its initial value (0), same rationale
        // as RodResiduals::computeResiduals.
        for (int k = m_neq_th + m_neq_irr; k < m_neq_j; ++k) rj[k] = yj[k];
    }
}

// ---------------------------------------------------------------------------
// computeGlobalUpdate — direct-substitution solve of the 3 global rows
// (fggen backward-Euler ODE, fgrel AE, gpres AE). See header: none of these
// three equations implicitly depends on its own unknown given the current
// layer states, so no Newton solve is needed here, unlike per-layer blocks.
// Mirrors computeGlobalResiduals' RHS formulas exactly, solved directly.
// ---------------------------------------------------------------------------
void FredMNaResiduals::updateSubchannelField(double t) const {
    if (m_subchannel == nullptr) return;
    std::vector<double> qqv_all(m_geom.nz);
    for (int j = 0; j < m_geom.nz; ++j)
        qqv_all[j] = layerPower(j, t);
    m_subchannel->updateCoolantField(t, qqv_all.data(), m_geom.dz0.data());
}

void FredMNaResiduals::computeGlobalUpdate(double t, double dt, double fggen_old,
                                            double& fggen_new, double& fgrel_new,
                                            double& gpres_new) const {
    m_Tplenum = m_plenumTFn ? m_plenumTFn(t) : layerCoolant(m_geom.nz - 1, t);

    double fgGrate = 0.0;
    for (int j = 0; j < m_geom.nz; ++j) {
        const auto& s = m_state.layers[j];
        const double vol       = m_geom.azf0 * m_geom.dz0[j];
        const double rate_fiss = s.qqv / (200.0 * 1.60214e-13);
        fgGrate += 0.25 * rate_fiss * vol / MNA_AVOGADRO;
    }
    fgGrate = std::max(fgGrate, 0.0);

    fggen_new = fggen_old + dt * fgGrate;
    fgrel_new = m_fgrel_grsis;
    gpres_new = computeGasPressure(m_state);
}

} // namespace fred
