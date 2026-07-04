#include "FredMNaSolver.hpp"
#include "FredMNaIrradiationPhysics.hpp"
#include "FredMNaGrsis.hpp"
#include "FredMNaGapBehavior.hpp"
#include "FredMNaFailure.hpp"
#include "FredMNaDenseSolve.hpp"
#include "platform/Constants.hpp"
#include "platform/RodResiduals.hpp"

#include "hdf5.h"
// No SUNDIALS IDA integrator here (see class comment in FredMNaSolver.hpp) —
// only nvector_serial (plain N_Vector data container, reused for
// checkpoint/HDF5-restart compatibility with FredSolverBase) and
// sundials_context (required to allocate an N_Vector at all).
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>

#include <cmath>
#include <stdexcept>
#include <memory>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cstring>
#include <cassert>
#include <chrono>

namespace fred {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
FredMNaSolver::FredMNaSolver(const FuelRodGeometry& geom,
                               UPuZr&                 fuel,
                               HT9&                   clad)
    : FredSolverBase(geom),
      m_fuel  (fuel),
      m_clad  (clad),
      m_gap_na(),
      m_res   (geom, fuel, clad, m_gap_na, 0.0),
      m_state (geom.nf, geom.nc, geom.nz),
      m_powerTab  (geom.nz),
      m_coolantTab(geom.nz)
{}

// ---------------------------------------------------------------------------
// Initial gas pressure setter
// ---------------------------------------------------------------------------
void FredMNaSolver::setInitialGasPressure(double p_MPa) {
    m_res.setInitialGasPressure(p_MPa);
}

// ---------------------------------------------------------------------------
// Boundary-condition hooks
// ---------------------------------------------------------------------------
void FredMNaSolver::onPowerDensityPerLayerSet(const std::vector<TimeTable>& tabs) {
    for (int j = 0; j < m_geom.nz; ++j)
        if (j < (int)tabs.size()) m_powerTab[j] = tabs[j];
}

void FredMNaSolver::onCoolantTemperaturePerLayerSet(const std::vector<TimeTable>&) {
    throw std::logic_error(
        "FredMNaSolver uses subchannel calculation for coolant temperature and pressure\n"
        "with a Robin boundary condition for the cladding outer surface heat flux.\n\n"
        "Example: solver.setCoolantChannel(dhyd, xarea, flowr, T_inlet_times, T_inlet_vals)\n\n"
        "See FredMNaSubchannelMode.cpp for implementation details.");
}

void FredMNaSolver::setCoolantChannel(
        double dhyd, double xarea, double flowr,
        std::vector<double> T_inlet_times,
        std::vector<double> T_inlet_vals,
        FredMNaSubchannelMode::HtcCorrelation corr)
{
    TimeTable inlet_tab(std::move(T_inlet_times), std::move(T_inlet_vals));
    auto inlet_fn = [inlet_tab](double t) { return inlet_tab(t); };
    m_subchannel_mode = std::make_unique<FredMNaSubchannelMode>(
        dhyd, xarea, flowr, std::move(inlet_fn),
        m_geom.rfo0, m_na_props, m_geom.nz, corr);
    m_res.setSubchannelMode(*m_subchannel_mode);
}

void FredMNaSolver::onCoolantPressureSet(const TimeTable& tbl) {
    m_res.setCoolantPressure([tbl](double t) { return tbl(t); });
}

void FredMNaSolver::onInitialTemperatureSet(double T0) {
    m_res.setInitialTemperature(T0);
}

void FredMNaSolver::onPlenumTemperatureSet(const TimeTable& tbl) {
    m_plenumTab = tbl;
}

// ---------------------------------------------------------------------------
// initState
// ---------------------------------------------------------------------------
void FredMNaSolver::initState() {
    const int nf = m_geom.nf, nc = m_geom.nc, nz = m_geom.nz;
    const double pu = m_fuel.puContent();
    const double zr = m_fuel.zrContent();

    // Register per-layer time functions with the residuals.
    {
        std::vector<TimeSeries> powerFns(nz), coolantFns(nz);
        for (int j = 0; j < nz; ++j) {
            auto ptab = m_powerTab[j];
            powerFns[j]   = [ptab](double t) { return ptab(t); };
            auto ctab = m_coolantTab[j];
            coolantFns[j] = [ctab](double t) { return ctab(t); };
        }
        m_res.setLayerPowerFunctions  (std::move(powerFns));
        m_res.setLayerCoolantFunctions(std::move(coolantFns));
    }
    if (!m_plenumTab.empty()) {
        auto ptab = m_plenumTab;
        m_res.setPlenumTemperature([ptab](double t) { return ptab(t); });
    }

    m_state.fggen        = 0.0;
    m_state.fgrel        = 0.0;
    m_state.bupave_FIMA  = 0.0;

    for (int j = 0; j < nz; ++j) {
        auto& s = m_state.layers[j];
        std::fill(s.T.begin(),    s.T.end(),    m_T0);
        std::fill(s.dT.begin(),   s.dT.end(),   0.0);
        s.qqv  = m_powerTab[j](0.0);
        s.bup  = 0.0; s.bup_FIMA = 0.0; s.buhard_FIMA = 0.0;
        s.gap  = m_geom.rci0 - m_geom.rfo0;
        s.pfc  = 0.0; s.hgap = 0.0;
        s.gapOpen = s.gap > (m_geom.ruff + m_geom.rufc);
        s.flag = s.gapOpen ? "open" : "clos";
        std::fill(s.efs.begin(),  s.efs.end(),  0.0);
        std::fill(s.defs.begin(), s.defs.end(), 0.0);

        for (int i = 0; i < nf; ++i) {
            auto& nd = s.nodes[i];
            nd.zr_wf = zr; nd.pu_wf = pu; nd.ur_wf = 1.0 - pu - zr;
            nd.phase = "alpha"; nd.pfrac = 1.0;
            nd.psod  = 0.0; nd.poros_tot = 0.0; nd.poros_gas = 0.0;

            const double dens0 = m_fuel.referenceDensity();
            const double r0    = m_geom.rad0[i];
            const double r1    = m_geom.rad0[i+1];
            nd.dvol  = PI * (r1*r1 - r0*r0) * m_geom.dz0[j];
            nd.mass  = dens0 * nd.dvol;
            const double ma_v = (zr / 91.22 + pu / 244.06 + (1.0-pu-zr) / 238.02891);
            const double ma   = 1.0 / ma_v;
            nd.zr_at = (zr / 91.22) * ma;
            nd.pu_at = (pu / 244.06) * ma;
            nd.ur_at = ((1.0-pu-zr) / 238.02891) * ma;
            const double tot_atoms = nd.mass * MNA_AVOGADRO / (ma * 1.0e-3);
            nd.c_zr  = tot_atoms * nd.zr_at / nd.dvol;
        }

        s.xwast = 0.0; s.clfuel = 0.0;
        const double ma = 1.0 / ((1.0-pu-zr) / 238.02891e-3 + pu / 244.06e-3 + zr / 91.22e-3);
        s.ntot = m_fuel.referenceDensity() * MNA_AVOGADRO / ma;

        std::fill(s.eft.begin(),  s.eft.end(),  0.0);
        std::fill(s.efh.begin(),  s.efh.end(),  0.0);
        std::fill(s.efr.begin(),  s.efr.end(),  0.0);
        s.efz = 0.0;
        std::fill(s.sigfh.begin(), s.sigfh.end(), 0.0);
        std::fill(s.sigfr.begin(), s.sigfr.end(), 0.0);
        std::fill(s.sigfz.begin(), s.sigfz.end(), 0.0);
        std::fill(s.et.begin(),   s.et.end(),   0.0);
        std::fill(s.eh.begin(),   s.eh.end(),   0.0);
        std::fill(s.er.begin(),   s.er.end(),   0.0);
        s.ez = 0.0;
        std::fill(s.sigh.begin(), s.sigh.end(), 0.0);
        std::fill(s.sigr.begin(), s.sigr.end(), 0.0);
        std::fill(s.sigz.begin(), s.sigz.end(), 0.0);
        for (int i = 0; i < nf+nc; ++i) s.rad[i] = m_geom.rad0[i];
        s.rfi = m_geom.rfi0; s.rfo = m_geom.rfo0;
        s.rci = m_geom.rci0; s.rco = m_geom.rco0;
    }

    m_res.initAlgebraicState(m_state);
}

// ---------------------------------------------------------------------------
// storeOutput
// ---------------------------------------------------------------------------
void FredMNaSolver::storeOutput(double t) {
    const int nf = m_geom.nf, nc = m_geom.nc, nz = m_geom.nz;

    m_t_prev = t;

    m_times_out.push_back(t);
    m_gpres_out.push_back(m_state.gpres);
    m_fggen_out.push_back(m_state.fggen);
    m_fgrel_out.push_back(m_state.fgrel);
    for (int j = 0; j < nz; ++j)
        for (int i = 0; i < nf+nc; ++i)
            m_T_out.push_back(m_state.layers[j].T[i]);
    double gap_avg = 0.0, bup_avg = 0.0;
    for (int j = 0; j < nz; ++j) {
        gap_avg += m_state.layers[j].gap;
        bup_avg += m_state.layers[j].bup;
    }
    m_gap_out.push_back(gap_avg / nz);
    m_bup_out.push_back(bup_avg / nz);

    // Per-layer coolant T, gap, and hgap at this output step
    for (int j = 0; j < nz; ++j) {
        m_T_co_out.push_back(m_subchannel_mode ? m_subchannel_mode->T_co(j) : 0.0);
        m_gap_layer_out.push_back(m_state.layers[j].gap);
        m_hgap_layer_out.push_back(m_state.layers[j].hgap);
    }

    double xwmax = 0.0;
    for (int j = 0; j < nz; ++j)
        xwmax = std::max(xwmax, m_state.layers[j].xwast);
    m_xwast_out.push_back(xwmax);

    // GRSIS swelling (spatially averaged over layers and nodes)
    double swtot_avg = 0.0, swopen_avg = 0.0;
    for (int j = 0; j < nz; ++j)
        for (int i = 0; i < nf-1; ++i) {
            swtot_avg  += m_state.layers[j].nodes[i].grsis.swtot;
            swopen_avg += m_state.layers[j].nodes[i].grsis.swopen;
        }
    const double n_nodes = std::max(nz * (nf - 1), 1);
    m_swtot_out .push_back(swtot_avg  / n_nodes);
    m_swopen_out.push_back(swopen_avg / n_nodes);

    // Failure criteria (worst layer)
    {
        const double pu = m_fuel.puContent();
        std::vector<double> T_fuel(nz), T_clav(nz), sig_eff(nz), zr_wf(nz);
        for (int j = 0; j < nz; ++j) {
            const auto& s = m_state.layers[j];
            T_fuel[j] = s.T[0];
            T_clav[j] = 0.5 * (s.T[nf] + s.T[nf+nc-1]);
            sig_eff[j] = 0.0;
            double area_sum = 0.0;
            for (int k = 0; k < nc; ++k) {
                const double da = s.rad[nf+k+1]*s.rad[nf+k+1] - s.rad[nf+k]*s.rad[nf+k];
                sig_eff[j] += s.sig[k] * da;
                area_sum   += da;
            }
            if (area_sum > 0.0) sig_eff[j] /= area_sum;
            sig_eff[j] *= 1.0e6; // MPa → Pa

            double zr_sum = 0.0;
            for (int i = 0; i < nf; ++i) zr_sum += s.nodes[i].zr_wf;
            zr_wf[j] = (nf > 0) ? zr_sum / nf : m_fuel.zrContent();
        }
        FredMNaFailureState fcrit;
        computeFailureCriteria(m_clad, m_fuel, nz,
                               T_fuel.data(), T_clav.data(),
                               sig_eff.data(), zr_wf.data(), fcrit);

        double burst_max = 0.0, melt_max = 0.0;
        for (int j = 0; j < nz; ++j) {
            burst_max = std::max(burst_max, fcrit.fcrit_burst[j]);
            melt_max  = std::max(melt_max,  fcrit.fcrit_melt[j]);
        }
        m_burst_out.push_back(burst_max);
        m_melt_out .push_back(melt_max);
    }

    storeYYP();
    flushH5Step(t);
}

// ---------------------------------------------------------------------------
// HDF5 streaming — MNA-specific
// ---------------------------------------------------------------------------
void FredMNaSolver::openAppH5Datasets()
{
    hid_t file  = (hid_t)m_h5ctx.file;
    hid_t therm = (hid_t)m_h5ctx.grp_therm;
    auto mk1d = [](hid_t loc, const char* name) -> hid_t {
        hsize_t dims[1]={0}, maxd[1]={H5S_UNLIMITED}, ck[1]={64};
        hid_t sp=H5Screate_simple(1,dims,maxd);
        hid_t pr=H5Pcreate(H5P_DATASET_CREATE); H5Pset_chunk(pr,1,ck);
        hid_t ds=H5Dcreate(loc,name,H5T_NATIVE_DOUBLE,sp,H5P_DEFAULT,pr,H5P_DEFAULT);
        H5Pclose(pr); H5Sclose(sp); return ds;
    };
    m_h5mna.ds_gap = mk1d(therm, "gap_width");
    hid_t bup = H5Gcreate(file, "burnup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_h5mna.grp_burnup = bup;
    m_h5mna.ds_gpres   = mk1d(bup, "gpres");
    m_h5mna.ds_fggen   = mk1d(bup, "fggen");
    m_h5mna.ds_fgrel   = mk1d(bup, "fgrel");
    m_h5mna.ds_bup     = mk1d(bup, "bup");
    m_h5mna.ds_xwast   = mk1d(bup, "xwast");
    m_h5mna.ds_swtot   = mk1d(bup, "swtot");
    m_h5mna.ds_swopen  = mk1d(bup, "swopen");
    m_h5mna.ds_burst   = mk1d(bup, "burst");
    m_h5mna.ds_melt    = mk1d(bup, "melt");
}

void FredMNaSolver::appendAppH5Row()
{
    if (m_h5mna.ds_gap < 0) return;
    auto app = [](hid_t ds, hsize_t step, double val) {
        hsize_t n=step+1; H5Dset_extent(ds,&n);
        hid_t fsp=H5Dget_space(ds); hsize_t st=step,ct=1;
        H5Sselect_hyperslab(fsp,H5S_SELECT_SET,&st,nullptr,&ct,nullptr);
        hid_t msp=H5Screate_simple(1,&ct,nullptr);
        H5Dwrite(ds,H5T_NATIVE_DOUBLE,msp,fsp,H5P_DEFAULT,&val);
        H5Sclose(msp); H5Sclose(fsp);
    };
    const hsize_t s = (hsize_t)m_h5ctx.nsteps;
    app((hid_t)m_h5mna.ds_gap,    s, m_gap_out   .empty() ? 0.0 : m_gap_out   .back());
    app((hid_t)m_h5mna.ds_gpres,  s, m_gpres_out .empty() ? 0.0 : m_gpres_out .back());
    app((hid_t)m_h5mna.ds_fggen,  s, m_fggen_out .empty() ? 0.0 : m_fggen_out .back());
    app((hid_t)m_h5mna.ds_fgrel,  s, m_fgrel_out .empty() ? 0.0 : m_fgrel_out .back());
    app((hid_t)m_h5mna.ds_bup,    s, m_bup_out   .empty() ? 0.0 : m_bup_out   .back());
    app((hid_t)m_h5mna.ds_xwast,  s, m_xwast_out .empty() ? 0.0 : m_xwast_out .back());
    app((hid_t)m_h5mna.ds_swtot,  s, m_swtot_out .empty() ? 0.0 : m_swtot_out .back());
    app((hid_t)m_h5mna.ds_swopen, s, m_swopen_out.empty() ? 0.0 : m_swopen_out.back());
    app((hid_t)m_h5mna.ds_burst,  s, m_burst_out .empty() ? 0.0 : m_burst_out .back());
    app((hid_t)m_h5mna.ds_melt,   s, m_melt_out  .empty() ? 0.0 : m_melt_out  .back());
}

void FredMNaSolver::trimAppOutputVectors()
{
    auto trim1 = [](std::vector<double>& v){ if(v.size()>1){ double x=v.back(); v={x}; } };
    trim1(m_gpres_out); trim1(m_fggen_out); trim1(m_fgrel_out);
    trim1(m_gap_out);   trim1(m_bup_out);   trim1(m_xwast_out);
    trim1(m_swtot_out); trim1(m_swopen_out);
    trim1(m_burst_out); trim1(m_melt_out);
}

void FredMNaSolver::closeAppH5Datasets()
{
    if (m_h5mna.ds_gap    >= 0) H5Dclose((hid_t)m_h5mna.ds_gap);
    if (m_h5mna.ds_gpres  >= 0) H5Dclose((hid_t)m_h5mna.ds_gpres);
    if (m_h5mna.ds_fggen  >= 0) H5Dclose((hid_t)m_h5mna.ds_fggen);
    if (m_h5mna.ds_fgrel  >= 0) H5Dclose((hid_t)m_h5mna.ds_fgrel);
    if (m_h5mna.ds_bup    >= 0) H5Dclose((hid_t)m_h5mna.ds_bup);
    if (m_h5mna.ds_xwast  >= 0) H5Dclose((hid_t)m_h5mna.ds_xwast);
    if (m_h5mna.ds_swtot  >= 0) H5Dclose((hid_t)m_h5mna.ds_swtot);
    if (m_h5mna.ds_swopen >= 0) H5Dclose((hid_t)m_h5mna.ds_swopen);
    if (m_h5mna.ds_burst  >= 0) H5Dclose((hid_t)m_h5mna.ds_burst);
    if (m_h5mna.ds_melt   >= 0) H5Dclose((hid_t)m_h5mna.ds_melt);
    if (m_h5mna.grp_burnup>= 0) H5Gclose((hid_t)m_h5mna.grp_burnup);
    m_h5mna = MnaH5Ctx{};
}

// ---------------------------------------------------------------------------
// logStepOutput
// ---------------------------------------------------------------------------
void FredMNaSolver::logStepOutput(double tret, double dt_next) {
    std::cout << "  t = " << tret << " s   next dt = " << dt_next
              << " s   gpres = " << m_state.gpres
              << " MPa  fgrel = " << m_state.fgrel * R_GAS * 293.15 / 1e5 * 1e6
              << " cm3  bupave = " << m_state.bupave_FIMA * 100 << " at%\n";
}

// ---------------------------------------------------------------------------
// runTimeLoop hooks
// ---------------------------------------------------------------------------
void FredMNaSolver::unpackCurrentState() {
    m_res.unpack(N_VGetArrayPointer((N_Vector)m_y),
                 N_VGetArrayPointer((N_Vector)m_yp), m_state);
}

void FredMNaSolver::handleGapClosed(int j) {
    FredSolverBase::handleGapClosed(j);
    m_state.layers[j].buhard_FIMA = m_state.layers[j].bup_FIMA;
    m_state.layers[j].flag = "clos";
}

void FredMNaSolver::handleGapReopened(int j) {
    // Unreachable post gap-behaviour-model refactor (Decision D4): no IDA
    // root-finder is registered for FRED-M-Na any more (setupGapRoots is not
    // called in run()), so runTimeLoop's root-event dispatch — the only
    // caller of this override — never fires for this app. Kept as a no-op
    // guard rather than deleted outright: setting flag="open" here would
    // violate the monotonic open->soft->clos ratchet if this were ever
    // (incorrectly) reached.
    assert(false && "FredMNaSolver::handleGapReopened is unreachable: "
                     "no gap root-finder is registered for FRED-M-Na");
}

// ---------------------------------------------------------------------------
// afterAcceptedStep -- per-step physics: burnup, GRSIS, Zr redistribution,
//                      cladding wastage, fgrel from GRSIS cgrel.
// Called by runOneStepLoop after every accepted backward-Euler step (this
// app does not use SUNDIALS IDA -- see class comment above).
// ---------------------------------------------------------------------------
void FredMNaSolver::afterAcceptedStep(double t, double dt) {
    if (dt <= 0.0) return;

    const int nf = m_geom.nf, nz = m_geom.nz;

    // Snapshot volumetric swelling (efs) from the previous accepted step,
    // before unpackCurrentState() overwrites m_state with the newly accepted
    // backward-Euler state. FredMNaGapBehavior::update needs both old and new swtot
    // (= 3*efs) per node to compute this step's ΔSwtot for the directional
    // swelling partition (Baseir.for anisotropy, Decision D8).
    std::vector<std::vector<double>> efs_prev(nz, std::vector<double>(nf));
    for (int j = 0; j < nz; ++j)
        for (int i = 0; i < nf; ++i)
            efs_prev[j][i] = m_state.layers[j].efs[i];

    unpackCurrentState();  // sync m_state from the accepted backward-Euler y/yp

    const double pu = m_fuel.puContent();
    const double zr = m_fuel.zrContent();
    const GrsisParams grsis_p = (m_grsis_mode == GrsisDataMode::FEAST)
                                ? GrsisParams::feast() : GrsisParams::grsis();

    const double molwe = 1.0 / (pu / 0.244 + (1.0 - zr - pu) / 0.238029);
    const double hm    = m_fuel.referenceDensity() * MNA_AVOGADRO / molwe;

    // Burnup must be updated for all layers before the anisotropy factor
    // below, which picks the current peak-linear-heat-rate layer.
    for (int j = 0; j < nz; ++j)
        m_state.layers[j].bup_FIMA += m_state.layers[j].qqv / (3.204354e-11 * hm) * dt;

    // Fuel axial anisotropy factor (rod-scalar, upuzrFanis / Baseir.for): drives
    // early localized "soft" fuel-clad contact from anisotropic metal-fuel
    // swelling, well before the bulk gap geometrically closes. Recomputed from
    // the current peak-power layer, frozen once that layer's burnup exceeds
    // 0.5% FIMA (matches legacy).
    {
        int anis_max = 0;
        double ql_max = -1.0;
        const double area = PI * (m_geom.rfo0*m_geom.rfo0 - m_geom.rfi0*m_geom.rfi0);
        for (int j = 0; j < nz; ++j) {
            const double ql = m_state.layers[j].qqv * area;
            if (ql > ql_max) { ql_max = ql; anis_max = j; }
        }
        if (m_state.layers[anis_max].bup_FIMA < 0.5e-2) {
            m_state.fanis_F    = ql_max / (m_geom.rfo0 * 2.0) * 1.0e-4;
            m_state.fanis_coef = upuzrFanis(m_state.fanis_F, pu);
        }
    }
    const double gap0  = m_geom.rci0 - m_geom.rfo0;
    const double rough = m_geom.ruff + m_geom.rufc;

    for (int j = 0; j < nz; ++j) {
        auto& s = m_state.layers[j];

        // Gap-open geometry check (Decision D3): FRED-M-Na does not use a
        // SUNDIALS root-finder for gap events (it doesn't use SUNDIALS IDA
        // at all any more, see class comment) — replicate legacy's step-
        // based detection by checking gap geometry directly at each fixed
        // step. m_res.setGapOpen keeps the residuals' GapStateManager (which
        // drives both the mechanics gap-contact bookkeeping and the thermal
        // gap-conductance model) in sync with this solver-side state.
        //
        // Unlike the old IDA-based scheme, a discontinuous BC change here
        // needs no special re-initialisation step: the one-step integrator
        // has no adaptive-step/multistep history to invalidate — the next
        // fixed backward-Euler step's Newton solve just sees the new
        // equations directly, exactly like every other step (and exactly
        // like legacy, which also has no such notion).
        //
        // Once closed, this must never re-evaluate: Baseir.for (line 546-549)
        // `goto`s straight past its entire gap/flag/pfc block once
        // flag=='clos', so pfc stays pinned at the last contact value
        // forever for the rest of the run (its reopening check at line 570 is
        // explicitly gated `fmat.ne.'upuzr'`, i.e. never taken for metal
        // fuel). Recomputing `s.gap > rough` unconditionally every step (as
        // this used to) instead re-opens/re-closes the gap dozens of times
        // once `s.gap` settles to within float noise of `rough` post-closure
        // (observed from ~t=1330d on layers 21-23), each flip spuriously
        // toggling the mechanical contact BC and gap conductance — this was
        // the direct cause of the gas-pressure discontinuity around t=1500d.
        if (s.flag != "clos") {
            const bool new_gapOpen = (s.gap > rough);
            if (new_gapOpen != s.gapOpen)
                m_res.setGapOpen(j, new_gapOpen);
            s.gapOpen = new_gapOpen;
        }

        // Gap-contact ratchet ("open" -> "soft" -> "clos") + directional
        // swelling partition (Baseir.for), owned by FredMNaGapBehavior
        // (Decision D7): monotonic — legacy never reopens a soft/closed gap
        // for upuzr fuel (Baseir.for's reopening check is gated on
        // fmat.ne.'upuzr').
        {
            std::vector<double> swtot_new(nf), swtot_old(nf);
            for (int i = 0; i < nf; ++i) {
                swtot_new[i] = 3.0 * s.efs[i];
                swtot_old[i] = 3.0 * efs_prev[j][i];
            }
            FredMNaGapBehavior::update(s, gap0, rough, m_state.fanis_coef,
                                        swtot_new.data(), swtot_old.data(), nf);
        }

        for (int i = 0; i < nf; ++i) {
            auto& nd = s.nodes[i];
            auto ph = upuzrPhase(s.T[i], pu, zr);
            nd.phase = ph.phase; nd.pfrac = ph.pfrac;
            nd.psod  = sodiumInfiltration(s.bup_FIMA, s.buhard_FIMA,
                                           m_geom.rad0[i], m_geom.rad0[nf-1],
                                           nd.phase, s.flag);
            nd.poros_tot = std::max(0.0, s.efs[i]);
            nd.poros_gas = 0.5 * nd.poros_tot;
        }

        // Irradiation conductivity correction
        {
            double T_avg = 0.0;
            for (int i = 0; i < nf; ++i) T_avg += s.T[i];
            T_avg /= nf;
            const auto& nd0 = s.nodes[std::min(nf/2, nf-1)];
            double k_irr = m_fuel.thermalConductivityIrradiated(
                T_avg, s.bup_FIMA,
                nd0.poros_tot, nd0.poros_gas, nd0.psod);
            double k_fresh = m_fuel.thermalConductivity(T_avg);
            s.k_irr_factor = (k_fresh > 1e-6) ? k_irr / k_fresh : 1.0;
        }

        if (m_enable_zr && nf > 1) {
            std::vector<double> T_ctr(nf-1);
            for (int i = 0; i < nf-1; ++i)
                T_ctr[i] = 0.5 * (s.T[i] + s.T[i+1]);

            std::vector<std::string> phase_ctr(nf-1);
            std::vector<double> pfrac_ctr(nf-1);
            for (int i = 0; i < nf-1; ++i) {
                auto ph = upuzrPhase(T_ctr[i], pu, zr);
                phase_ctr[i] = ph.phase;
                pfrac_ctr[i] = ph.pfrac;
            }

            std::vector<double> zr_at(nf-1), pu_at(nf-1), ur_at(nf-1);
            std::vector<double> zr_wf(nf-1), pu_wf(nf-1), ur_wf(nf-1);
            std::vector<double> zr_atoms(nf-1), pu_atoms(nf-1), ur_atoms(nf-1);
            std::vector<double> c_zr_arr(nf-1), mass_arr(nf-1), dvol_arr(nf-1);

            for (int i = 0; i < nf-1; ++i) {
                zr_at[i]    = 0.5*(s.nodes[i].zr_at + s.nodes[i+1].zr_at);
                pu_at[i]    = 0.5*(s.nodes[i].pu_at + s.nodes[i+1].pu_at);
                ur_at[i]    = 0.5*(s.nodes[i].ur_at + s.nodes[i+1].ur_at);
                zr_wf[i]    = 0.5*(s.nodes[i].zr_wf + s.nodes[i+1].zr_wf);
                pu_wf[i]    = 0.5*(s.nodes[i].pu_wf + s.nodes[i+1].pu_wf);
                ur_wf[i]    = 0.5*(s.nodes[i].ur_wf + s.nodes[i+1].ur_wf);
                c_zr_arr[i] = 0.5*(s.nodes[i].c_zr  + s.nodes[i+1].c_zr);
                dvol_arr[i] = s.nodes[i].dvol;
                mass_arr[i] = s.nodes[i].mass;
                const double tot    = zr_at[i] + pu_at[i] + ur_at[i];
                const double ma_v   = (tot > 0) ? (zr_at[i]*91.22 + pu_at[i]*244.06 + ur_at[i]*238.02891) : 91.22;
                zr_atoms[i] = c_zr_arr[i] * dvol_arr[i];
                pu_atoms[i] = pu_at[i] / std::max(tot, 1e-30) * mass_arr[i] * MNA_AVOGADRO / (ma_v * 1e-3);
                ur_atoms[i] = ur_at[i] / std::max(tot, 1e-30) * mass_arr[i] * MNA_AVOGADRO / (ma_v * 1e-3);
            }

            upuzrZirconiumRedistribution(
                nf-1, T_ctr.data(), m_geom.rad0.data(),
                phase_ctr.data(), pfrac_ctr.data(),
                zr_at.data(), pu_at.data(), ur_at.data(),
                zr_wf.data(), pu_wf.data(), ur_wf.data(),
                zr_atoms.data(), pu_atoms.data(), ur_atoms.data(),
                c_zr_arr.data(), mass_arr.data(), dvol_arr.data(), dt);

            for (int i = 0; i < nf-1; ++i) {
                s.nodes[i].zr_at = zr_at[i];
                s.nodes[i].pu_at = pu_at[i];
                s.nodes[i].ur_at = ur_at[i];
                s.nodes[i].zr_wf = zr_wf[i];
                s.nodes[i].pu_wf = pu_wf[i];
                s.nodes[i].ur_wf = ur_wf[i];
                s.nodes[i].c_zr  = c_zr_arr[i];
            }
        }

        // GRSIS bubble swelling model (per radial node)
        if (m_enable_grsis && nf > 1) {
            const double fissd = s.qqv / (200.0 * 1.60214e-13);
            const double bup0  = s.bup_FIMA - s.qqv / (3.204354e-11 * hm) * dt;
            const double gpres_Pa = m_state.gpres * 1.0e6;

            // upuzrGrsis is a 1:1 port of the legacy Fortran GRSIS_metalfuel,
            // which integrates the bubble/gas-release ODEs with unguarded
            // explicit Euler and no internal substepping. The legacy driver
            // only ever calls it with dt = dtout (~1000 s in the reference
            // input deck). Here dt is FredMNaSolver's own fixed outer step
            // (up to dtout, which can be set to days/months by the caller) —
            // far outside the regime the explicit-Euler scheme is stable in
            // (cgrel diverges past cggen) if applied directly at that size.
            // Substep internally to stay near the legacy cadence, same
            // pattern as ht9CladWastage below.
            constexpr double GRSIS_MAX_SUBSTEP = 1000.0; // s
            const int    n_sub  = std::max(1, (int)std::ceil(dt / GRSIS_MAX_SUBSTEP));
            const double sub_dt = dt / n_sub;

            for (int i = 0; i < nf-1; ++i) {
                auto& nd = s.nodes[i];
                const double pfc_Pa = s.pfc * 1.0e6;
                const double swtot_before = nd.grsis.swtot;
                for (int k = 0; k < n_sub; ++k) {
                    const double bup_k0 = bup0 + (s.bup_FIMA - bup0) * double(k)   / n_sub;
                    const double bup_k1 = bup0 + (s.bup_FIMA - bup0) * double(k+1) / n_sub;
                    upuzrGrsis(nd.grsis, m_grsis_first, sub_dt,
                               s.T[i], fissd, bup_k1, bup_k0,
                               gpres_Pa, pfc_Pa,
                               s.sigfr[i], s.sigfh[i], s.sigfz[i],
                               !s.gapOpen, grsis_p);
                }

                nd.poros_tot = nd.grsis.swtot;
                nd.poros_gas = nd.grsis.frtot;

                // Couple GRSIS's total (solid + gas-bubble) swelling into the
                // mechanical efs ODE: d(efs)/dt = d(swtot)/dt / 3 (volumetric ->
                // linear strain, same convention as the empirical swsol-only rate
                // it replaces). One-step lag, same explicit-coupling pattern as
                // fgrel (setFgrelFromGrsis) above.
                if (dt > 0.0)
                    s.defs_grsis[i] = (nd.grsis.swtot - swtot_before) / (3.0 * dt);
            }
            // Outermost node mirrors the last inner node
            s.nodes[nf-1].poros_tot = s.nodes[nf-2].poros_tot;
            s.nodes[nf-1].poros_gas = s.nodes[nf-2].poros_gas;
            s.nodes[nf-1].grsis     = s.nodes[nf-2].grsis;
            s.defs_grsis[nf-1]      = s.defs_grsis[nf-2];
        }

        if (m_enable_waste) {
            std::vector<double> qqv_arr  (1, s.qqv);
            std::vector<double> cit_arr  (1, s.T[nf]);
            std::vector<std::string> flag_arr(1, s.flag);
            std::vector<double> ntot_arr (1, s.ntot);
            std::vector<double> clf_arr  (1, s.clfuel);
            std::vector<double> xw_arr   (1, s.xwast);
            std::vector<double> dz0_arr  (1, m_geom.dz0[j]);

            ht9CladWastage(1, dt, std::max(t, 1.0),
                            qqv_arr.data(), cit_arr.data(), flag_arr.data(),
                            m_fuel.referenceDensity(), pu, zr,
                            m_geom.rfo0, dz0_arr.data(),
                            ntot_arr.data(), clf_arr.data(), xw_arr.data());
            s.clfuel = clf_arr[0];
            s.xwast  = xw_arr[0];
            s.ntot   = ntot_arr[0];
        }

        // Push this step's irradiation physics (Zr redistribution, GRSIS,
        // k_irr_factor, defs_grsis, flag, efsz/efsh/efsr, ...) into the
        // residuals' own state — see FredMNaResiduals::syncAuxLayerState for
        // why this is needed.
        m_res.syncAuxLayerState(j, s);
    }

    m_state.bupave_FIMA = 0.0;
    for (int j = 0; j < nz; ++j) m_state.bupave_FIMA += m_state.layers[j].bup_FIMA;
    m_state.bupave_FIMA /= nz;

    // Compute fgrel from GRSIS cgrel (AE: replaces empirical ODE).
    if (m_enable_grsis && nf > 1) {
        double fgrel_new = 0.0;
        for (int j = 0; j < nz; ++j)
            for (int i = 0; i < nf-1; ++i)
                fgrel_new += m_state.layers[j].nodes[i].grsis.cgrel
                           * m_state.layers[j].nodes[i].dvol;
        fgrel_new /= MNA_AVOGADRO;
        m_state.fgrel = std::max(0.0, fgrel_new);
        m_res.setFgrelFromGrsis(m_state.fgrel);
    }

    m_grsis_first = false;
}

// ---------------------------------------------------------------------------
// run
// ---------------------------------------------------------------------------
void FredMNaSolver::run(double tend, double dtout, bool all_steps) {
    m_times_out.clear(); m_T_out.clear();
    m_gpres_out.clear(); m_fggen_out.clear(); m_fgrel_out.clear();
    m_gap_out.clear();   m_bup_out.clear();
    m_xwast_out.clear(); m_y_out.clear(); m_yp_out.clear();
    m_swtot_out.clear(); m_swopen_out.clear();
    m_burst_out.clear(); m_melt_out.clear();
    m_T_co_out.clear(); m_gap_layer_out.clear(); m_hgap_layer_out.clear();
    openH5File("fred-m-na");

    const bool restarting    = (m_restart_time >= 0.0);
    const bool is_snapshot   = m_restart_is_snapshot;
    const double t_start     = restarting ? m_restart_time : 0.0;
    m_restart_time        = -1.0;   // consumed
    m_restart_is_snapshot = false;  // consumed

    if (!restarting) {
        m_t_prev      = 0.0;
        m_grsis_first = true;
        initState();
    } else {
        // m_t_prev and m_grsis_first were restored by readAppCheckpoint.
        // Re-register per-layer profiles with the residuals assembler.
        const int nz = m_geom.nz;
        std::vector<TimeSeries> powerFns(nz), coolantFns(nz);
        for (int j = 0; j < nz; ++j) {
            auto ptab = m_powerTab[j];
            powerFns[j]   = [ptab](double t) { return ptab(t); };
            auto ctab = m_coolantTab[j];
            coolantFns[j] = [ctab](double t) { return ctab(t); };
        }
        m_res.setLayerPowerFunctions  (std::move(powerFns));
        m_res.setLayerCoolantFunctions(std::move(coolantFns));
        if (!m_plenumTab.empty()) {
            auto ptab = m_plenumTab;
            m_res.setPlenumTemperature([ptab](double t) { return ptab(t); });
        }
        // Restore elapsed time in residuals.
        // For snapshot restarts: elapsed time was already restored by readAppCheckpoint;
        // do not override it with t_start=0 which would reset accumulated irradiation history.
        if (!is_snapshot)
            m_res.setElapsedTime(t_start);
    }

    std::cout << "FRED-M-Na: neq=" << m_res.neqTotal() << ", nz=" << m_geom.nz
              << ", nf=" << m_geom.nf << ", nc=" << m_geom.nc << "\n";

    // ---- One-step integrator setup: N_Vector as a plain data container only
    // (no SUNDIALS IDA integrator — see class comment in FredMNaSolver.hpp).
    // Kept as N_Vector (not raw std::vector<double>) purely so the existing
    // checkpoint/snapshot/HDF5-restart plumbing in FredSolverBase — which
    // reads m_y/m_yp directly — keeps working unmodified. ----
    const int neq = m_res.neqTotal();
    SUNContext sunctx;
    if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != 0)
        throw std::runtime_error("FredMNaSolver: SUNContext_Create failed");
    m_sunctx = sunctx;
    N_Vector y  = N_VNew_Serial(neq, sunctx);
    N_Vector yp = N_VNew_Serial(neq, sunctx);
    if (!y || !yp) throw std::runtime_error("FredMNaSolver: N_VNew_Serial failed");
    m_y  = y;
    m_yp = yp;
    std::fill(N_VGetArrayPointer(yp), N_VGetArrayPointer(yp) + neq, 0.0);

    if (!restarting) {
        // No IDACalcIC-equivalent "consistent IC" pass: legacy has none
        // either (Baseir.for starts its Picard loop directly from the
        // initialised state), and the first fixed backward-Euler step below
        // resolves the full nonlinear system (algebraic + differential rows
        // together) exactly like every subsequent step, so any small
        // inconsistency in initState()/initAlgebraicState()'s seed values is
        // absorbed there rather than needing a separate phase.
        m_res.pack(m_state, N_VGetArrayPointer(y));
    } else {
        applyRestartOneStep(t_start);
    }

    // Seed FredMNaResiduals' own internal m_state (globals + all layers) once
    // up front: computeLayerResiduals only ever updates ONE layer per call
    // (that's the point — O(1), not O(nz)), so without this the very first
    // outer Picard sweep's computeGlobalUpdate call would see zero-
    // initialized layer geometry (see setGlobalState/primeState doc comments).
    m_res.primeState(N_VGetArrayPointer(y), N_VGetArrayPointer(yp));

    m_t_prev = t_start;
    unpackCurrentState();
    storeOutput(t_start);

    const double dt_step = (m_step_size > 0.0) ? m_step_size : dtout;
    std::cout << "FRED-M-Na: integrating to t=" << tend << " s"
              << " (one-step backward-Euler, fixed dt=" << dt_step << " s)"
              << (all_steps ? " (all-steps mode)" : "")
              << (restarting ? " (restarting from t=" + std::to_string(t_start) + " s)" : "")
              << "\n";

    runOneStepLoop(tend, dtout, all_steps, t_start);
    closeH5File();
    std::cout << "FRED-M-Na: done\n";
}

// ---------------------------------------------------------------------------
// applyRestartOneStep — inject checkpoint state without touching IDA
// (FredSolverBase::applyRestartToIDA calls IDAReInit, which needs m_ida_mem —
// never created in one-step mode). Mirrors applyRestartToIDA minus that call.
// ---------------------------------------------------------------------------
void FredMNaSolver::applyRestartOneStep(double /*t_start*/) {
    const int nz = m_geom.nz;
    double* y  = N_VGetArrayPointer((N_Vector)m_y);
    double* yp = N_VGetArrayPointer((N_Vector)m_yp);
    std::copy(m_restart_y .begin(), m_restart_y .end(), y);
    std::copy(m_restart_yp.begin(), m_restart_yp.end(), yp);
    for (int j = 0; j < nz; ++j)
        rodResiduals().restoreGapState(j, m_restart_gap_open[j], m_restart_axial_offset[j]);
}

// ---------------------------------------------------------------------------
// runOneStepLoop — fixed-step, backward-Euler, "always accept" time loop.
// Replaces FredSolverBase::runTimeLoop for this app (see class comment).
// Mirrors runTimeLoop's output/checkpoint/snapshot bookkeeping, but drives
// the step with solveStepBackwardEuler() instead of IDASolve(IDA_ONE_STEP),
// and has no root-finding / gap-event dispatch: gap-state transitions are
// detected and applied directly inside afterAcceptedStep every fixed step
// (see there), with no re-initialisation step needed afterwards.
// ---------------------------------------------------------------------------
void FredMNaSolver::runOneStepLoop(double tend, double dtout, bool all_steps, double t_start) {
    double t = t_start;
    double t_next_out = t_start + dtout;
    const double dt_req = (m_step_size > 0.0) ? m_step_size : dtout;

    std::vector<double> pending_snaps = m_snapshot_timings;
    m_snapshot_count = 0;

    while (t < tend - 1.0e-12 * dtout) {
        double dt = std::min(dt_req, t_next_out - t);
        dt = std::min(dt, tend - t);
        if (dt <= 0.0) break;
        const double t_new = t + dt;

        solveStepBackwardEuler(t_new, dt);

        const auto t_phys_start = std::chrono::steady_clock::now();
        afterAcceptedStep(t_new, dt);
        if (m_verbosity >= 5) {
            const double phys_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t_phys_start).count();
            std::cout << "      [physics] afterAcceptedStep wall=" << phys_ms << " ms\n";
        }
        t = t_new;

        const bool at_dtout = (t >= t_next_out - 1.0e-12 * dtout);
        bool saved_snap_this_step = false;

        if (at_dtout) {
            unpackCurrentState();
            storeOutput(t);
            if (m_verbosity >= 3)
                logStepOutput(t, dt_req);

            if (!m_ckpt_prefix.empty())
                saveCheckpoint(t);

            if (!m_snapshot_prefix.empty()) {
                for (auto it = pending_snaps.begin(); it != pending_snaps.end(); ) {
                    if (std::abs(t - *it) < 0.5 * dtout) {
                        ++m_snapshot_count;
                        std::string fname = m_snapshot_prefix + "_frame"
                            + std::to_string(m_snapshot_count) + ".snapshot";
                        saveSnapshot(t, fname);
                        saved_snap_this_step = true;
                        it = pending_snaps.erase(it);
                    } else {
                        ++it;
                    }
                }
            }

            t_next_out += dtout;
        } else if (all_steps) {
            unpackCurrentState();
            storeOutput(t);
        }

        if (t >= tend - 1.0e-12 * dtout) {
            if (!m_snapshot_prefix.empty() && !saved_snap_this_step) {
                ++m_snapshot_count;
                std::string fname = m_snapshot_prefix + "_frame"
                    + std::to_string(m_snapshot_count) + ".snapshot";
                saveSnapshot(t, fname);
            }
            break;
        }
    }
}

// ---------------------------------------------------------------------------
// solveStepBackwardEuler — advance y from t_new-dt to t_new via a fixed-dt
// backward-Euler outer Picard sweep: the 3 global rows (fggen/fgrel/gpres)
// are solved by direct substitution given the current layer states (none of
// them implicitly depends on its own unknown — see computeGlobalUpdate),
// then each axial layer's own local block is solved by dense Newton (layers
// only couple through the globals, never directly to each other — see
// computeLayerResiduals). The sweep is capped at MAX_OUTER iterations and
// ALWAYS accepted whatever it ends on — no rejection, no dt shrinking —
// replicating legacy FRED-M's "backward-Euler solution always accepted"
// Picard loop (Baseir.for: niter capped at 1000, 'goto 100' escape-and-
// continue on non-convergence) with a per-layer Newton solve in place of
// legacy's node-by-node successive substitution.
// ---------------------------------------------------------------------------
void FredMNaSolver::solveStepBackwardEuler(double t_new, double dt) {
    const int neq = m_res.neqTotal();
    double* y  = N_VGetArrayPointer((N_Vector)m_y);
    double* yp = N_VGetArrayPointer((N_Vector)m_yp);

    std::vector<double> y_old(y, y + neq);
    std::vector<double> id(neq);
    m_res.fillIdArray(id.data());

    const int nz    = m_geom.nz;
    const int neq_j = m_res.neqPerLayer();
    const double fggen_old = y_old[0];

    constexpr int    MAX_OUTER  = 60;
    constexpr double OUTER_RTOL = 1.0e-6;

    double gpres_prev = y[2];
    int    outer_used = 0;
    double layer_delta_last = 0.0;

    for (int outer = 0; outer < MAX_OUTER; ++outer) {
        ++outer_used;

        // fggen/fgrel are genuinely explicit (direct substitution; neither
        // implicitly depends on its own unknown) — solved here directly.
        // gpres is NOT: its defining equation has real cross-coupling to
        // each layer's own mechanics (fuel/clad radii, axial strain), so it
        // is solved JOINTLY with each layer below instead (newtonSolveLayer)
        // for tighter within-step coupling than a separate direct-
        // substitution pass would give. computeGlobalUpdate's own gpres
        // estimate is discarded; the call is still needed for its fggen/
        // fgrel outputs and its m_Tplenum side effect.
        // Refresh subchannel coolant field (T_co/HTC per layer) for t_new —
        // computeLayerResiduals is deliberately scoped to one layer and
        // does not do this itself (see updateSubchannelField doc comment).
        m_res.updateSubchannelField(t_new);

        double fggen_new, fgrel_new, gpres_hint;
        m_res.computeGlobalUpdate(t_new, dt, fggen_old, fggen_new, fgrel_new, gpres_hint);
        (void)gpres_hint;
        y[0] = fggen_new; y[1] = fgrel_new;
        yp[0] = (fggen_new - fggen_old) / dt; yp[1] = 0.0; yp[2] = 0.0;
        m_res.setGlobalState(fggen_new, fgrel_new, y[2]);

        double layer_delta = 0.0;
        for (int j = 0; j < nz; ++j) {
            double*       yj     = y  + 3 + j * neq_j;
            double*       ypj    = yp + 3 + j * neq_j;
            const double* yj_old = y_old.data() + 3 + j * neq_j;
            const double* idj    = id.data()    + 3 + j * neq_j;
            const double delta = newtonSolveLayer(j, t_new, dt, yj_old, idj, yj, ypj, neq_j,
                                                   fggen_new, fgrel_new, &y[2]);
            layer_delta = std::max(layer_delta, delta);
        }
        layer_delta_last = layer_delta;

        const double gpres_scale = std::max(std::fabs(y[2]), 1.0e-8);
        const double outer_err   = std::fabs(y[2] - gpres_prev) / gpres_scale;
        gpres_prev = y[2];

        if (outer_err < OUTER_RTOL && layer_delta < 1.0e-6)
            break;
        // Otherwise keep sweeping; if this was the last allowed outer
        // iteration, the loop simply ends and whatever y currently holds is
        // accepted as this step's result (no retry with a smaller dt).
    }

    if (m_verbosity >= 4) {
        std::cout << "    [onestep] t=" << t_new << " dt=" << dt
                  << " outer_iters=" << outer_used
                  << " layer_delta=" << layer_delta_last << "\n";
    }
}

// ---------------------------------------------------------------------------
// newtonSolveLayer — dense Newton solve of one axial layer's local block
// (backward-Euler discretisation: yp[k] = (yj[k]-yj_old[k])/dt for
// differential rows, 0 for algebraic rows, matching id[]/fillIdArray's
// convention). See FredMNaDenseSolve.hpp for the generic solver.
// ---------------------------------------------------------------------------
double FredMNaSolver::newtonSolveLayer(int j, double t_new, double dt,
                                        const double* yj_old, const double* idj,
                                        double* yj, double* ypj, int n,
                                        double fggen_new, double fgrel_new,
                                        double* gpres_inout) const {
    // +1 unknown: gpres, solved jointly with this layer's own mechanics
    // (see header doc comment — gpres's own defining equation
    // gpres = computeGasPressure(...) has real dependence on this layer's
    // fuel/clad radii and axial strain unknowns, so solving it together
    // gives tighter within-step coupling).
    const int m = n + 1;
    std::vector<double> y_ext(m);
    std::copy(yj, yj + n, y_ext.begin());
    y_ext[n] = *gpres_inout;

    std::vector<double> fd_scale(m);
    m_res.fillFdScale(fd_scale.data());
    fd_scale[n] = 1.0; // gpres [MPa]

    auto residualFn = [&](const double* ytrial, double* rout) {
        std::vector<double> ypl(n);
        for (int k = 0; k < n; ++k)
            ypl[k] = (idj[k] > 0.5) ? (ytrial[k] - yj_old[k]) / dt : 0.0;
        m_res.setGlobalState(fggen_new, fgrel_new, ytrial[n]);
        m_res.computeLayerResiduals(j, t_new, ytrial, ypl.data(), rout);
        rout[n] = ytrial[n] - m_res.computeGasPressureCurrent();
    };

    const double delta = denseNewtonSolve(m, residualFn, y_ext.data(), /*max_iter=*/25,
                                           /*atol=*/1.0e-7, fd_scale.data());

    std::copy(y_ext.begin(), y_ext.begin() + n, yj);
    *gpres_inout = y_ext[n];
    for (int k = 0; k < n; ++k)
        ypj[k] = (idj[k] > 0.5) ? (yj[k] - yj_old[k]) / dt : 0.0;
    m_res.setGlobalState(fggen_new, fgrel_new, *gpres_inout);

    return delta;
}

// ---------------------------------------------------------------------------
// writeAppCheckpoint — FRED-M-Na specific state beyond y/yp/gap.
// ---------------------------------------------------------------------------
void FredMNaSolver::writeAppCheckpoint(std::ostream& os) const {
    const int nf = m_geom.nf, nz = m_geom.nz;

    // Elapsed irradiation time and m_t_prev
    double elapsed = m_res.elapsedTime();
    os.write(reinterpret_cast<const char*>(&elapsed),  8);
    os.write(reinterpret_cast<const char*>(&m_t_prev), 8);
    const double fgrel_grsis = m_state.fgrel;
    os.write(reinterpret_cast<const char*>(&fgrel_grsis), 8);
    os.write(reinterpret_cast<const char*>(&m_state.fanis_F),    8);
    os.write(reinterpret_cast<const char*>(&m_state.fanis_coef), 8);

    // App checkpoint format version — local to FredMNaSolver's own payload,
    // independent of the shared "FREDCKPT" header version in FredSolverBase
    // (which stays untouched, Decision D1). v2 (gap-behaviour-model refactor,
    // Decision D6): added per-layer flag + directional swelling
    // (efsz/efsh/efsr); breaking change from v1, acceptable per D6 (no
    // long-running M-Na jobs depend on old checkpoint files).
    const int32_t mna_ckpt_ver = 2;
    os.write(reinterpret_cast<const char*>(&mna_ckpt_ver), 4);

    // Per-layer state
    for (int j = 0; j < nz; ++j) {
        const auto& s = m_state.layers[j];

        // 5 scalars per layer
        os.write(reinterpret_cast<const char*>(&s.bup_FIMA),    8);
        os.write(reinterpret_cast<const char*>(&s.buhard_FIMA), 8);
        os.write(reinterpret_cast<const char*>(&s.xwast),       8);
        os.write(reinterpret_cast<const char*>(&s.clfuel),      8);
        os.write(reinterpret_cast<const char*>(&s.ntot),        8);

        // Gap contact flag ("open"/"soft"/"clos"), fixed 8 chars (padded with nulls)
        {
            char flag_buf[8] = {};
            std::strncpy(flag_buf, s.flag.c_str(), 7);
            os.write(flag_buf, 8);
        }

        // Directional swelling accumulators (Baseir.for anisotropy), nf doubles each
        os.write(reinterpret_cast<const char*>(s.efsz.data()), nf * 8);
        os.write(reinterpret_cast<const char*>(s.efsh.data()), nf * 8);
        os.write(reinterpret_cast<const char*>(s.efsr.data()), nf * 8);

        // Per-node state
        for (int i = 0; i < nf; ++i) {
            const auto& nd = s.nodes[i];
            // Composition: 9 doubles
            os.write(reinterpret_cast<const char*>(&nd.zr_wf), 8);
            os.write(reinterpret_cast<const char*>(&nd.pu_wf), 8);
            os.write(reinterpret_cast<const char*>(&nd.ur_wf), 8);
            os.write(reinterpret_cast<const char*>(&nd.zr_at), 8);
            os.write(reinterpret_cast<const char*>(&nd.pu_at), 8);
            os.write(reinterpret_cast<const char*>(&nd.ur_at), 8);
            os.write(reinterpret_cast<const char*>(&nd.c_zr),  8);
            os.write(reinterpret_cast<const char*>(&nd.mass),  8);
            os.write(reinterpret_cast<const char*>(&nd.dvol),  8);
            // Phase/porosity: 4 doubles
            os.write(reinterpret_cast<const char*>(&nd.pfrac),     8);
            os.write(reinterpret_cast<const char*>(&nd.psod),      8);
            os.write(reinterpret_cast<const char*>(&nd.poros_tot), 8);
            os.write(reinterpret_cast<const char*>(&nd.poros_gas), 8);
            // Phase string: fixed 8 chars (padded with nulls)
            char phase_buf[8] = {};
            std::strncpy(phase_buf, nd.phase.c_str(), 7);
            os.write(phase_buf, 8);
            // GRSIS state: 18 doubles
            const auto& g = nd.grsis;
            os.write(reinterpret_cast<const char*>(&g.cgfm),    8);
            os.write(reinterpret_cast<const char*>(&g.cgb1),    8);
            os.write(reinterpret_cast<const char*>(&g.cgb2),    8);
            os.write(reinterpret_cast<const char*>(&g.cgb31),   8);
            os.write(reinterpret_cast<const char*>(&g.cgb32),   8);
            os.write(reinterpret_cast<const char*>(&g.cggen),   8);
            os.write(reinterpret_cast<const char*>(&g.cgrel),   8);
            os.write(reinterpret_cast<const char*>(&g.rb1),     8);
            os.write(reinterpret_cast<const char*>(&g.rb2),     8);
            os.write(reinterpret_cast<const char*>(&g.vb1),     8);
            os.write(reinterpret_cast<const char*>(&g.vb2),     8);
            os.write(reinterpret_cast<const char*>(&g.swclose), 8);
            os.write(reinterpret_cast<const char*>(&g.swopen),  8);
            os.write(reinterpret_cast<const char*>(&g.swgas),   8);
            os.write(reinterpret_cast<const char*>(&g.swsol),   8);
            os.write(reinterpret_cast<const char*>(&g.swtot),   8);
            os.write(reinterpret_cast<const char*>(&g.frg),     8);
            os.write(reinterpret_cast<const char*>(&g.frtot),   8);
        }
    }
}

// ---------------------------------------------------------------------------
// readAppCheckpoint — restore FRED-M-Na specific state from checkpoint.
// ---------------------------------------------------------------------------
void FredMNaSolver::readAppCheckpoint(std::istream& is) {
    const int nf = m_geom.nf, nz = m_geom.nz;

    // Elapsed irradiation time and m_t_prev
    double elapsed;
    is.read(reinterpret_cast<char*>(&elapsed),  8);
    is.read(reinterpret_cast<char*>(&m_t_prev), 8);
    double fgrel_grsis;
    is.read(reinterpret_cast<char*>(&fgrel_grsis), 8);
    is.read(reinterpret_cast<char*>(&m_state.fanis_F),    8);
    is.read(reinterpret_cast<char*>(&m_state.fanis_coef), 8);
    m_res.setElapsedTime(elapsed);
    m_res.setFgrelFromGrsis(fgrel_grsis);
    m_grsis_first = false;  // GRSIS already initialised

    int32_t mna_ckpt_ver = 0;
    is.read(reinterpret_cast<char*>(&mna_ckpt_ver), 4);
    if (mna_ckpt_ver != 2)
        throw std::runtime_error(
            "FredMNaSolver::readAppCheckpoint: unsupported checkpoint version "
            "(expected 2, got " + std::to_string(mna_ckpt_ver) +
            ") — old-format M-Na checkpoints predate the gap-behaviour-model "
            "refactor and are not supported (Decision D6)");

    // Per-layer state (m_state already has correct size from constructor)
    for (int j = 0; j < nz; ++j) {
        auto& s = m_state.layers[j];

        is.read(reinterpret_cast<char*>(&s.bup_FIMA),    8);
        is.read(reinterpret_cast<char*>(&s.buhard_FIMA), 8);
        is.read(reinterpret_cast<char*>(&s.xwast),       8);
        is.read(reinterpret_cast<char*>(&s.clfuel),      8);
        is.read(reinterpret_cast<char*>(&s.ntot),        8);

        {
            char flag_buf[8] = {};
            is.read(flag_buf, 8);
            s.flag = std::string(flag_buf);
        }
        is.read(reinterpret_cast<char*>(s.efsz.data()), nf * 8);
        is.read(reinterpret_cast<char*>(s.efsh.data()), nf * 8);
        is.read(reinterpret_cast<char*>(s.efsr.data()), nf * 8);

        for (int i = 0; i < nf; ++i) {
            auto& nd = s.nodes[i];
            is.read(reinterpret_cast<char*>(&nd.zr_wf), 8);
            is.read(reinterpret_cast<char*>(&nd.pu_wf), 8);
            is.read(reinterpret_cast<char*>(&nd.ur_wf), 8);
            is.read(reinterpret_cast<char*>(&nd.zr_at), 8);
            is.read(reinterpret_cast<char*>(&nd.pu_at), 8);
            is.read(reinterpret_cast<char*>(&nd.ur_at), 8);
            is.read(reinterpret_cast<char*>(&nd.c_zr),  8);
            is.read(reinterpret_cast<char*>(&nd.mass),  8);
            is.read(reinterpret_cast<char*>(&nd.dvol),  8);
            is.read(reinterpret_cast<char*>(&nd.pfrac),     8);
            is.read(reinterpret_cast<char*>(&nd.psod),      8);
            is.read(reinterpret_cast<char*>(&nd.poros_tot), 8);
            is.read(reinterpret_cast<char*>(&nd.poros_gas), 8);
            char phase_buf[8] = {};
            is.read(phase_buf, 8);
            nd.phase = std::string(phase_buf);
            auto& g = nd.grsis;
            is.read(reinterpret_cast<char*>(&g.cgfm),    8);
            is.read(reinterpret_cast<char*>(&g.cgb1),    8);
            is.read(reinterpret_cast<char*>(&g.cgb2),    8);
            is.read(reinterpret_cast<char*>(&g.cgb31),   8);
            is.read(reinterpret_cast<char*>(&g.cgb32),   8);
            is.read(reinterpret_cast<char*>(&g.cggen),   8);
            is.read(reinterpret_cast<char*>(&g.cgrel),   8);
            is.read(reinterpret_cast<char*>(&g.rb1),     8);
            is.read(reinterpret_cast<char*>(&g.rb2),     8);
            is.read(reinterpret_cast<char*>(&g.vb1),     8);
            is.read(reinterpret_cast<char*>(&g.vb2),     8);
            is.read(reinterpret_cast<char*>(&g.swclose), 8);
            is.read(reinterpret_cast<char*>(&g.swopen),  8);
            is.read(reinterpret_cast<char*>(&g.swgas),   8);
            is.read(reinterpret_cast<char*>(&g.swsol),   8);
            is.read(reinterpret_cast<char*>(&g.swtot),   8);
            is.read(reinterpret_cast<char*>(&g.frg),     8);
            is.read(reinterpret_cast<char*>(&g.frtot),   8);
        }
    }
}

} // namespace fred
