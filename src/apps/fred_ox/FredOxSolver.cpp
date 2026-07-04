#include "FredOxSolver.hpp"
#include "fuelpelletmaterial/FredOxMOX.hpp"
#include "platform/CladdingMaterial.hpp"
#include "platform/Constants.hpp"

#include "hdf5.h"
#include <ida/ida.h>
#include <nvector/nvector_serial.h>

#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>

namespace fred {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
FredOxSolver::FredOxSolver(const FuelRodGeometry& geom,
                             FredOxMOX&             mox,
                             const CladdingMaterial& clad,
                             FredOxGapMaterial&     gap)
    : FredIdaSolverBase(geom),
      m_mox  (mox),
      m_clad (clad),
      m_gap  (gap),
      m_res  (geom, mox, clad, gap,
              0.0,
              mox.referenceDensity(),
              mox.puContent(),
              mox.stoichiometry(),
              1.0),
      m_state(geom.nf, geom.nc, geom.nz),
      m_powerTab  (geom.nz),
      m_coolantTab(geom.nz)
{}

// ---------------------------------------------------------------------------
// Boundary-condition hooks
// ---------------------------------------------------------------------------
void FredOxSolver::onPowerDensityPerLayerSet(const std::vector<TimeTable>& tabs) {
    for (int j = 0; j < m_geom.nz; ++j)
        if (j < (int)tabs.size()) m_powerTab[j] = tabs[j];
}

void FredOxSolver::onCoolantTemperaturePerLayerSet(const std::vector<TimeTable>& tabs) {
    for (int j = 0; j < m_geom.nz; ++j)
        if (j < (int)tabs.size()) m_coolantTab[j] = tabs[j];
}

void FredOxSolver::onCoolantPressureSet(const TimeTable& tbl) {
    m_res.setCoolantPressure([tbl](double t) { return tbl(t); });
}

void FredOxSolver::onInitialTemperatureSet(double T0) {
    m_res.setInitialTemperature(T0);
}

void FredOxSolver::onPlenumTemperatureSet(const TimeTable& tbl) {
    m_plenumTab = tbl;
}

// ---------------------------------------------------------------------------
// OX-only boundary condition setters
// ---------------------------------------------------------------------------
void FredOxSolver::setInitialGasPressure(double gpres0_MPa) {
    m_res.setInitialGasPressure(gpres0_MPa);
    const double mu0 = gpres0_MPa * 1.0e6 * m_res.vfree0() / T_REF / R_GAS;
    m_res.setGasInventory(mu0);
    m_gap.setGasInventory(mu0);
    m_state.mu0 = mu0;
}

void FredOxSolver::setSwellingMultiplier(double fswelmlt) {
    m_res.setSwellingMultiplier(fswelmlt);
}

// ---------------------------------------------------------------------------
// initState
// ---------------------------------------------------------------------------
void FredOxSolver::initState() {
    const int nf = m_geom.nf, nc = m_geom.nc, nz = m_geom.nz;

    // Register per-layer functions with the residuals.
    {
        std::vector<TimeSeries> powerFns(nz), coolantFns(nz);
        for (int j = 0; j < nz; ++j) {
            auto ptab = m_powerTab[j];
            powerFns[j] = [ptab](double t) { return ptab(t); };
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

    m_state.fggen = 0.0;
    m_state.fgrel = 0.0;

    for (int j = 0; j < nz; ++j) {
        auto& s = m_state.layers[j];
        std::fill(s.T.begin(),     s.T.end(),     m_T0);
        std::fill(s.dT.begin(),    s.dT.end(),    0.0);
        s.qqv = m_powerTab[j](0.0);
        s.bup = 0.0;
        s.gap = m_geom.rci0 - m_geom.rfo0;
        s.pfc = 0.0;  s.hgap = 0.0;
        s.gapOpen = s.gap > (m_geom.ruff + m_geom.rufc);
        std::fill(s.efs.begin(),   s.efs.end(),   0.0);
        std::fill(s.eft.begin(),   s.eft.end(),   0.0);
        std::fill(s.efh.begin(),   s.efh.end(),   0.0);
        std::fill(s.efr.begin(),   s.efr.end(),   0.0);
        s.efz = 0.0;
        std::fill(s.sigfh.begin(), s.sigfh.end(), 0.0);
        std::fill(s.sigfr.begin(), s.sigfr.end(), 0.0);
        std::fill(s.sigfz.begin(), s.sigfz.end(), 0.0);
        std::fill(s.et.begin(),    s.et.end(),    0.0);
        std::fill(s.eh.begin(),    s.eh.end(),    0.0);
        std::fill(s.er.begin(),    s.er.end(),    0.0);
        s.ez = 0.0;
        std::fill(s.sigh.begin(),  s.sigh.end(),  0.0);
        std::fill(s.sigr.begin(),  s.sigr.end(),  0.0);
        std::fill(s.sigz.begin(),  s.sigz.end(),  0.0);
        for (int i = 0; i < nf+nc; ++i) s.rad[i] = m_geom.rad0[i];
        s.rfi = m_geom.rfi0;  s.rfo = m_geom.rfo0;
        s.rci = m_geom.rci0;  s.rco = m_geom.rco0;
    }
    m_res.initAlgebraicState(m_state);
}

// ---------------------------------------------------------------------------
// storeOutput
// ---------------------------------------------------------------------------
void FredOxSolver::storeOutput(double t) {
    const int nf = m_geom.nf, nc = m_geom.nc, nz = m_geom.nz;
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
    storeYYP();
    flushH5Step(t);
}

// ---------------------------------------------------------------------------
// HDF5 streaming — OX-specific
// ---------------------------------------------------------------------------
void FredOxSolver::openAppH5Datasets()
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
    m_h5ox.ds_gap = mk1d(therm, "gap_width");
    hid_t bup = H5Gcreate(file, "burnup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_h5ox.grp_burnup = bup;
    m_h5ox.ds_gpres   = mk1d(bup, "gpres");
    m_h5ox.ds_fggen   = mk1d(bup, "fggen");
    m_h5ox.ds_fgrel   = mk1d(bup, "fgrel");
    m_h5ox.ds_bup     = mk1d(bup, "bup");
}

void FredOxSolver::appendAppH5Row()
{
    if (m_h5ox.ds_gap < 0) return;
    auto app = [](hid_t ds, hsize_t step, double val) {
        hsize_t n=step+1; H5Dset_extent(ds,&n);
        hid_t fsp=H5Dget_space(ds); hsize_t st=step,ct=1;
        H5Sselect_hyperslab(fsp,H5S_SELECT_SET,&st,nullptr,&ct,nullptr);
        hid_t msp=H5Screate_simple(1,&ct,nullptr);
        H5Dwrite(ds,H5T_NATIVE_DOUBLE,msp,fsp,H5P_DEFAULT,&val);
        H5Sclose(msp); H5Sclose(fsp);
    };
    const hsize_t s = (hsize_t)m_h5ctx.nsteps;
    app((hid_t)m_h5ox.ds_gap,   s, m_gap_out  .empty() ? 0.0 : m_gap_out  .back());
    app((hid_t)m_h5ox.ds_gpres, s, m_gpres_out.empty() ? 0.0 : m_gpres_out.back());
    app((hid_t)m_h5ox.ds_fggen, s, m_fggen_out.empty() ? 0.0 : m_fggen_out.back());
    app((hid_t)m_h5ox.ds_fgrel, s, m_fgrel_out.empty() ? 0.0 : m_fgrel_out.back());
    app((hid_t)m_h5ox.ds_bup,   s, m_bup_out  .empty() ? 0.0 : m_bup_out  .back());
}

void FredOxSolver::trimAppOutputVectors()
{
    auto trim1 = [](std::vector<double>& v){ if(v.size()>1){ double x=v.back(); v={x}; } };
    trim1(m_gpres_out); trim1(m_fggen_out); trim1(m_fgrel_out);
    trim1(m_gap_out);   trim1(m_bup_out);
}

void FredOxSolver::closeAppH5Datasets()
{
    if (m_h5ox.ds_gap   >= 0) H5Dclose((hid_t)m_h5ox.ds_gap);
    if (m_h5ox.ds_gpres >= 0) H5Dclose((hid_t)m_h5ox.ds_gpres);
    if (m_h5ox.ds_fggen >= 0) H5Dclose((hid_t)m_h5ox.ds_fggen);
    if (m_h5ox.ds_fgrel >= 0) H5Dclose((hid_t)m_h5ox.ds_fgrel);
    if (m_h5ox.ds_bup   >= 0) H5Dclose((hid_t)m_h5ox.ds_bup);
    if (m_h5ox.grp_burnup>=0) H5Gclose((hid_t)m_h5ox.grp_burnup);
    m_h5ox = OxH5Ctx{};
}

// ---------------------------------------------------------------------------
// runTimeLoop hooks
// ---------------------------------------------------------------------------
void FredOxSolver::unpackCurrentState() {
    m_res.unpack(N_VGetArrayPointer((N_Vector)m_y),
                 N_VGetArrayPointer((N_Vector)m_yp), m_state);
}

void FredOxSolver::logStepOutput(double tret, double dt_next) {
    std::cout << "  t = " << tret << " s   next dt = " << dt_next
              << " s   gpres = " << m_state.gpres
              << " MPa  fgrel = " << m_state.fgrel * R_GAS * 293.15 / 1e5 * 1e6 << " cm3\n";
}

// ---------------------------------------------------------------------------
// run
// ---------------------------------------------------------------------------
void FredOxSolver::run(double tend, double dtout, bool all_steps, int threads) {
    setThreads(threads);
    m_res.setNumThreads(this->threads());
    m_times_out.clear(); m_T_out.clear();
    m_gpres_out.clear(); m_fggen_out.clear(); m_fgrel_out.clear();
    m_gap_out.clear();   m_bup_out.clear();
    m_y_out.clear();     m_yp_out.clear();
    openH5File("fred-ox");

    const bool restarting = (m_restart_time >= 0.0);
    const double t_start  = restarting ? m_restart_time : 0.0;
    m_restart_time = -1.0;  // consumed

    if (!restarting) {
        initState();
    } else {
        // Re-register per-layer profiles with residuals (needed after restart).
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
    }

    const int neq = m_res.neqTotal();
    std::cout << "FRED-OX: neq=" << neq << ", nz=" << m_geom.nz
              << ", nf=" << m_geom.nf << ", nc=" << m_geom.nc << "\n";

    setupSundials(20000, t_start);
    setupGapRoots(m_geom.nz, RodResiduals::gapRoot);

    if (!restarting) {
        double t_ic = std::min(dtout, tend);
        if (t_ic <= 0.0) t_ic = 1.0;
        calcIC(t_ic);
    } else {
        applyRestartToIDA(t_start);
    }

    if (m_hot_start && !restarting) runHotStart();

    unpackCurrentState();
    storeOutput(t_start);

    std::cout << "FRED-OX: starting time integration to t=" << tend << " s"
              << (all_steps ? " (all-steps mode)" : "")
              << (restarting ? " (restarting from t=" + std::to_string(t_start) + " s)" : "")
              << "\n";

    runTimeLoop(tend, dtout, all_steps, t_start);
    closeH5File();
    std::cout << "FRED-OX: done\n";
}

} // namespace fred
