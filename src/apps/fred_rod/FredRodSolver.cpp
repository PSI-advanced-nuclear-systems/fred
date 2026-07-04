#include "FredRodSolver.hpp"
#include "platform/FuelPelletMaterial.hpp"
#include "platform/CladdingMaterial.hpp"

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
FredRodSolver::FredRodSolver(const FuelRodGeometry&    geom,
                              const FuelPelletMaterial& fuel,
                              const CladdingMaterial&   clad,
                              const GapMaterial&        gap)
    : FredIdaSolverBase(geom),
      m_fuel(fuel),
      m_clad(clad),
      m_gap(gap),
      m_pressure_model(std::make_unique<GapPressureModel>(0.1, geom, gap)),
      m_res  (geom, fuel, clad, gap, 0.0),
      m_state(geom.nf, geom.nc, geom.nz)
{
    m_res.setGapPressureModel(*m_pressure_model);
}

// ---------------------------------------------------------------------------
// Public setters
// ---------------------------------------------------------------------------
void FredRodSolver::setInitialGasPressure(double gpres0_MPa) {
    m_pressure_model->setInitialGasPressure(gpres0_MPa);
}

// ---------------------------------------------------------------------------
// Boundary-condition hooks
// ---------------------------------------------------------------------------
void FredRodSolver::onPowerDensityPerLayerSet(const std::vector<TimeTable>& tabs) {
    std::vector<TimeSeries> fns;
    fns.reserve(tabs.size());
    for (const auto& t : tabs)
        fns.push_back([t](double time) { return t(time); });
    m_res.setLayerPowerFunctions(std::move(fns));
}

void FredRodSolver::onCoolantTemperaturePerLayerSet(const std::vector<TimeTable>& tabs) {
    std::vector<TimeSeries> fns;
    fns.reserve(tabs.size());
    for (const auto& t : tabs)
        fns.push_back([t](double time) { return t(time); });
    m_res.setLayerCoolantFunctions(std::move(fns));
}

void FredRodSolver::onCoolantPressureSet(const TimeTable& tbl) {
    m_res.setCoolantPressure([tbl](double t) { return tbl(t); });
}

void FredRodSolver::onInitialTemperatureSet(double T0) {
    m_res.setInitialTemperature(T0);
    if (m_pressure_model) m_pressure_model->setInitialTemperature(T0);
}

void FredRodSolver::onPlenumTemperatureSet(const TimeTable& tbl) {
    m_pressure_model->setPlenumTemperature([tbl](double t) { return tbl(t); });
}

// ---------------------------------------------------------------------------
// ROD-only boundary condition setter
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// fillSolverIdArray — ROD layout: gpres AE + per-layer fuel T ODEs + inner-clad T ODEs
// ---------------------------------------------------------------------------
void FredRodSolver::fillSolverIdArray(double* id, int /*neq*/) const {
    // id[0] = 0.0 already (global gpres AE).
    const int neq_j = m_res.neqPerLayer();
    for (int j = 0; j < m_geom.nz; ++j) {
        int base = 1 + j * neq_j;
        for (int i = 0; i < m_geom.nf; ++i)         id[base + 4 + i] = 1.0;
        for (int i = 0; i < m_geom.nc - 1; ++i)     id[base + 4 + m_geom.nf + i] = 1.0;
    }
}

// ---------------------------------------------------------------------------
// initState
// ---------------------------------------------------------------------------
void FredRodSolver::initState() {
    const int nf = m_geom.nf, nc = m_geom.nc, nz = m_geom.nz;
    for (int j = 0; j < nz; ++j) {
        auto& s = m_state.layers[j];
        std::fill(s.T.begin(),    s.T.end(),    m_T0);
        std::fill(s.dT.begin(),   s.dT.end(),   0.0);
        s.qqv = 0.0;  s.gap = m_geom.rci0 - m_geom.rfo0;  s.pfc = 0.0;  s.hgap = 0.0;
        s.gapOpen = (s.gap > (m_geom.ruff + m_geom.rufc));
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
        for (int i = 0; i < nf-1;  ++i) s.drf[i] = m_geom.drf0;
        for (int i = 0; i < nc-1;  ++i) s.drc[i] = m_geom.drc0;
        s.rfi = m_geom.rfi0; s.rfo = m_geom.rfo0;
        s.rci = m_geom.rci0; s.rco = m_geom.rco0;
    }
    m_res.initAlgebraicState(m_state, m_fuel, m_clad);
}

// ---------------------------------------------------------------------------
// storeOutput
// ---------------------------------------------------------------------------
void FredRodSolver::storeOutput(double t) {
    const int nf = m_geom.nf, nc = m_geom.nc, nz = m_geom.nz;
    m_times_out.push_back(t);
    for (int j = 0; j < nz; ++j)
        for (int i = 0; i < nf+nc; ++i)
            m_T_out.push_back(m_state.layers[j].T[i]);

    double sigh_avg = 0.0, rco_avg = 0.0;
    double gap_avg = 0.0, pfc_avg = 0.0, rfo_avg = 0.0, rci_avg = 0.0;
    const double rough_gap = m_geom.ruff + m_geom.rufc;
    for (int j = 0; j < nz; ++j) {
        sigh_avg += m_state.layers[j].sigh[nc-1];
        rco_avg  += m_state.layers[j].rco;
        gap_avg  += m_state.layers[j].gapOpen
                 ? (m_state.layers[j].rci - m_state.layers[j].rfo) : rough_gap;
        pfc_avg  += m_state.layers[j].pfc;
        rfo_avg  += m_state.layers[j].rfo;
        rci_avg  += m_state.layers[j].rci;
    }
    m_sigh_outer_out.push_back(sigh_avg / nz);
    m_rco_out       .push_back(rco_avg  / nz);
    m_gap_out       .push_back(gap_avg  / nz);
    m_pfc_out       .push_back(pfc_avg  / nz);
    m_rfo_out       .push_back(rfo_avg  / nz);
    m_rci_out       .push_back(rci_avg  / nz);

    storeYYP();
    flushH5Step(t);
}

// ---------------------------------------------------------------------------
// HDF5 streaming — ROD-specific
// ---------------------------------------------------------------------------
void FredRodSolver::openAppH5Datasets()
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

    m_h5rod.ds_gap  = mk1d(therm, "gap_width");

    hid_t mech = H5Gcreate(file, "mechanical", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_h5rod.grp_mech = mech;
    m_h5rod.ds_sigh  = mk1d(mech, "sigh_outer");
    m_h5rod.ds_rco   = mk1d(mech, "rco");
    m_h5rod.ds_pfc   = mk1d(mech, "pfc");
    m_h5rod.ds_rfo   = mk1d(mech, "rfo");
    m_h5rod.ds_rci   = mk1d(mech, "rci");
}

void FredRodSolver::appendAppH5Row()
{
    if (m_h5rod.ds_gap < 0) return;
    auto app = [](hid_t ds, hsize_t step, double val) {
        hsize_t n = step+1; H5Dset_extent(ds,&n);
        hid_t fsp=H5Dget_space(ds); hsize_t st=step,ct=1;
        H5Sselect_hyperslab(fsp,H5S_SELECT_SET,&st,nullptr,&ct,nullptr);
        hid_t msp=H5Screate_simple(1,&ct,nullptr);
        H5Dwrite(ds,H5T_NATIVE_DOUBLE,msp,fsp,H5P_DEFAULT,&val);
        H5Sclose(msp); H5Sclose(fsp);
    };
    const hsize_t s = (hsize_t)m_h5ctx.nsteps;
    app((hid_t)m_h5rod.ds_gap,  s, m_gap_out .empty() ? 0.0 : m_gap_out .back());
    app((hid_t)m_h5rod.ds_sigh, s, m_sigh_outer_out.empty() ? 0.0 : m_sigh_outer_out.back());
    app((hid_t)m_h5rod.ds_rco,  s, m_rco_out .empty() ? 0.0 : m_rco_out .back());
    app((hid_t)m_h5rod.ds_pfc,  s, m_pfc_out .empty() ? 0.0 : m_pfc_out .back());
    app((hid_t)m_h5rod.ds_rfo,  s, m_rfo_out .empty() ? 0.0 : m_rfo_out .back());
    app((hid_t)m_h5rod.ds_rci,  s, m_rci_out .empty() ? 0.0 : m_rci_out .back());
}

void FredRodSolver::trimAppOutputVectors()
{
    auto trim1 = [](std::vector<double>& v){ if(v.size()>1){ double x=v.back(); v={x}; } };
    trim1(m_sigh_outer_out); trim1(m_rco_out); trim1(m_gap_out);
    trim1(m_pfc_out); trim1(m_rfo_out); trim1(m_rci_out);
}

void FredRodSolver::closeAppH5Datasets()
{
    if (m_h5rod.ds_gap  >= 0) H5Dclose((hid_t)m_h5rod.ds_gap);
    if (m_h5rod.ds_sigh >= 0) H5Dclose((hid_t)m_h5rod.ds_sigh);
    if (m_h5rod.ds_rco  >= 0) H5Dclose((hid_t)m_h5rod.ds_rco);
    if (m_h5rod.ds_pfc  >= 0) H5Dclose((hid_t)m_h5rod.ds_pfc);
    if (m_h5rod.ds_rfo  >= 0) H5Dclose((hid_t)m_h5rod.ds_rfo);
    if (m_h5rod.ds_rci  >= 0) H5Dclose((hid_t)m_h5rod.ds_rci);
    if (m_h5rod.grp_mech>= 0) H5Gclose((hid_t)m_h5rod.grp_mech);
    m_h5rod = RodH5Ctx{};
}

// ---------------------------------------------------------------------------
// runTimeLoop hooks
// ---------------------------------------------------------------------------
void FredRodSolver::unpackCurrentState() {
    m_res.unpack(N_VGetArrayPointer((N_Vector)m_y),
                 N_VGetArrayPointer((N_Vector)m_yp), m_state);
}

void FredRodSolver::afterGapEventsHandled(double tret, double dtout) {
    if (m_res.isMechOn()) {
        // Re-solve mechanical equilibrium in the new gap-constraint set.
        unpackCurrentState();
        m_res.solveMechanicalIC(m_state, tret, m_rtol, m_atol);
        m_res.pack(m_state, N_VGetArrayPointer((N_Vector)m_y));
    }
    reinitAfterEvent(tret, dtout);
}

// ---------------------------------------------------------------------------
// run
// ---------------------------------------------------------------------------
void FredRodSolver::run(double tend, double dtout, bool all_steps, int threads) {
    setThreads(threads);
    m_res.setNumThreads(this->threads());
    m_times_out.clear(); m_T_out.clear();
    m_sigh_outer_out.clear(); m_rco_out.clear();
    m_gap_out.clear();   m_pfc_out.clear();
    m_rfo_out.clear();   m_rci_out.clear();
    m_y_out.clear();     m_yp_out.clear();
    openH5File("fred-rod");

    const bool restarting = (m_restart_time >= 0.0);
    const double t_start  = restarting ? m_restart_time : 0.0;
    m_restart_time = -1.0;  // consumed; reset so next run() is cold unless re-loaded

    if (!restarting) {
        initState();
        if (m_res.isMechOn())
            m_res.solveMechanicalIC(m_state, 0.0, m_rtol, m_atol);
    }
    // For a restart, profiles were already registered by BC setters before loadCheckpoint.

    // Quasi-static path (heat off): no differential equations → Newton only.
    if (!restarting) {
        int nd = 0;
        if (m_res.isHeatOn()) nd = m_geom.nz * (m_geom.nf + (m_geom.nc - 1));
        if (nd == 0) {
            std::cout << "FRED-ROD: quasi-static (all-algebraic) — Newton at each output time\n";
            storeOutput(0.0);
            double tout = dtout;
            while (tout <= tend + 0.5 * dtout) {
                m_res.solveMechanicalIC(m_state, tout, m_rtol, m_atol);
                storeOutput(tout);
                tout += dtout;
            }
            std::cout << "FRED-ROD: done\n";
            closeH5File();
            return;
        }
    }

    // ---- IDA path ----
    setupSundials(10000, t_start);

    if (m_res.isMechOn()) {
        std::cout << "FRED-ROD: enabling root detection (nz=" << m_geom.nz << ")\n";
        setupGapRoots(m_geom.nz, RodResiduals::gapRoot);
    } else {
        std::cout << "FRED-ROD: root detection disabled (mechanics off)\n";
    }

    if (!restarting) {
        double t_ic = (dtout > 0.0) ? dtout : 1.0;
        calcIC(t_ic);
    } else {
        applyRestartToIDA(t_start);
    }

    if (m_hot_start && !restarting) runHotStart();

    unpackCurrentState();
    storeOutput(t_start);

    std::cout << "FRED-ROD: starting time integration"
              << (all_steps ? " (all-steps mode)" : "")
              << (restarting ? " (restarting from t=" + std::to_string(t_start) + " s)" : "")
              << "\n";

    runTimeLoop(tend, dtout, all_steps, t_start);
    closeH5File();
    std::cout << "FRED-ROD: done\n";
}

} // namespace fred
