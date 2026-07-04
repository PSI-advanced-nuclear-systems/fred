#include "FredSolverBase.hpp"
#include "RodResiduals.hpp"

#include "hdf5.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <cmath>

// ---------------------------------------------------------------------------
// HDF5 streaming helpers (anonymous namespace — internal to this TU)
// ---------------------------------------------------------------------------
namespace {

static hid_t h5_create_1d(hid_t loc, const char* name, hsize_t chunk = 64)
{
    hsize_t dims[1] = {0}, maxd[1] = {H5S_UNLIMITED}, ck[1] = {chunk};
    hid_t sp = H5Screate_simple(1, dims, maxd);
    hid_t pr = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(pr, 1, ck);
    hid_t ds = H5Dcreate(loc, name, H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, pr, H5P_DEFAULT);
    H5Pclose(pr); H5Sclose(sp);
    return ds;
}

static hid_t h5_create_2d(hid_t loc, const char* name, hsize_t ncols, hsize_t chunk_rows = 16)
{
    hsize_t dims[2] = {0, ncols}, maxd[2] = {H5S_UNLIMITED, ncols}, ck[2] = {chunk_rows, ncols};
    hid_t sp = H5Screate_simple(2, dims, maxd);
    hid_t pr = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(pr, 2, ck);
    hid_t ds = H5Dcreate(loc, name, H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, pr, H5P_DEFAULT);
    H5Pclose(pr); H5Sclose(sp);
    return ds;
}

static void h5_append_scalar(hid_t ds, hsize_t step, double val)
{
    hsize_t new_size = step + 1;
    H5Dset_extent(ds, &new_size);
    hid_t fsp = H5Dget_space(ds);
    hsize_t start = step, count = 1;
    H5Sselect_hyperslab(fsp, H5S_SELECT_SET, &start, nullptr, &count, nullptr);
    hid_t msp = H5Screate_simple(1, &count, nullptr);
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, msp, fsp, H5P_DEFAULT, &val);
    H5Sclose(msp); H5Sclose(fsp);
}

static void h5_append_row(hid_t ds, hsize_t step, const double* data, hsize_t ncols)
{
    hsize_t new_dims[2] = {step + 1, ncols};
    H5Dset_extent(ds, new_dims);
    hid_t fsp = H5Dget_space(ds);
    hsize_t start[2] = {step, 0}, count[2] = {1, ncols};
    H5Sselect_hyperslab(fsp, H5S_SELECT_SET, start, nullptr, count, nullptr);
    hid_t msp = H5Screate_simple(2, count, nullptr);
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, msp, fsp, H5P_DEFAULT, data);
    H5Sclose(msp); H5Sclose(fsp);
}

static void h5_attr_str(hid_t loc, const char* name, const char* val)
{
    hid_t tp = H5Tcopy(H5T_C_S1);
    H5Tset_size(tp, strlen(val) + 1);
    hid_t sp = H5Screate(H5S_SCALAR);
    hid_t at = H5Acreate(loc, name, tp, sp, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(at, tp, val);
    H5Aclose(at); H5Sclose(sp); H5Tclose(tp);
}

static void h5_attr_int(hid_t loc, const char* name, int val)
{
    hid_t sp = H5Screate(H5S_SCALAR);
    hid_t at = H5Acreate(loc, name, H5T_NATIVE_INT, sp, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(at, H5T_NATIVE_INT, &val);
    H5Aclose(at); H5Sclose(sp);
}

static void h5_attr_int_update(hid_t loc, const char* name, int val)
{
    if (H5Aexists(loc, name) > 0) H5Adelete(loc, name);
    h5_attr_int(loc, name, val);
}

} // anonymous namespace

namespace fred {

// ---------------------------------------------------------------------------
// Construction / Destruction
// ---------------------------------------------------------------------------
FredSolverBase::FredSolverBase(const FuelRodGeometry& geom)
    : m_geom(geom)
{}

FredSolverBase::~FredSolverBase() {
    if (m_y)       N_VDestroy((N_Vector)m_y);
    if (m_yp)      N_VDestroy((N_Vector)m_yp);
    if (m_sunctx)  SUNContext_Free((SUNContext*)&m_sunctx);
}

// ---------------------------------------------------------------------------
// Physics toggles
// ---------------------------------------------------------------------------
void FredSolverBase::setEnableHeatConduction(bool b) { rodResiduals().setEnableHeatConduction(b); }
void FredSolverBase::setEnableStressStrain  (bool b) { rodResiduals().setEnableStressStrain(b); }

// ---------------------------------------------------------------------------
// Shared boundary-condition setters (NVI)
// ---------------------------------------------------------------------------
void FredSolverBase::setPowerDensityHistoryPerLayer(
    std::vector<double> times,
    std::vector<std::vector<double>> qqv_per_layer)
{
    std::vector<TimeTable> tabs;
    tabs.reserve(qqv_per_layer.size());
    for (auto& qv : qqv_per_layer)
        tabs.emplace_back(times, qv);
    onPowerDensityPerLayerSet(tabs);
}

void FredSolverBase::setPowerDensityHistory(std::vector<double> times, std::vector<double> qv) {
    // Broadcast: build a single TimeTable, then replicate it for every layer.
    TimeTable tbl(times, qv);
    std::vector<TimeTable> tabs(m_geom.nz, tbl);
    onPowerDensityPerLayerSet(tabs);
}

void FredSolverBase::setCoolantTemperaturePerLayer(
    std::vector<double> times,
    std::vector<std::vector<double>> T_per_layer)
{
    std::vector<TimeTable> tabs;
    tabs.reserve(T_per_layer.size());
    for (auto& T : T_per_layer)
        tabs.emplace_back(times, T);
    onCoolantTemperaturePerLayerSet(tabs);
}

void FredSolverBase::setCoolantTemperature(std::vector<double> times, std::vector<double> T_K) {
    TimeTable tbl(times, T_K);
    std::vector<TimeTable> tabs(m_geom.nz, tbl);
    onCoolantTemperaturePerLayerSet(tabs);
}

void FredSolverBase::setCoolantPressureHistory(std::vector<double> times, std::vector<double> pcool_MPa) {
    TimeTable tbl(std::move(times), std::move(pcool_MPa));
    onCoolantPressureSet(tbl);
}

void FredSolverBase::setCoolantPressure(double pcool_MPa) {
    TimeTable tbl({0.0}, {pcool_MPa});
    onCoolantPressureSet(tbl);
}

void FredSolverBase::setPlenumTemperatureHistory(std::vector<double> times, std::vector<double> T_K) {
    TimeTable tbl(std::move(times), std::move(T_K));
    onPlenumTemperatureSet(tbl);
}

void FredSolverBase::setInitialTemperature(double T0_K) {
    m_T0 = T0_K;
    onInitialTemperatureSet(T0_K);
}

// ---------------------------------------------------------------------------
// peakFuelTemperature
// ---------------------------------------------------------------------------
std::vector<double> FredSolverBase::peakFuelTemperature() const {
    const int nf  = m_geom.nf;
    const int nc  = m_geom.nc;
    const int nz  = m_geom.nz;
    const int stride = nz * (nf + nc);

    std::vector<double> peaks;
    peaks.reserve(m_times_out.size());
    for (size_t t = 0; t < m_times_out.size(); ++t) {
        double peak = 0.0;
        for (int j = 0; j < nz; ++j)
            for (int i = 0; i < nf; ++i) {
                double T = m_T_out[t * stride + j*(nf+nc) + i];
                if (T > peak) peak = T;
            }
        peaks.push_back(peak);
    }
    return peaks;
}

// ---------------------------------------------------------------------------
// storeYYP — append current y/yp arrays to m_y_out / m_yp_out.
// ---------------------------------------------------------------------------
void FredSolverBase::storeYYP() {
    if (!m_y || !m_yp) return;
    const int neq = getSolverNeq();
    const double* ydata  = N_VGetArrayPointer((N_Vector)m_y);
    const double* ypdata = N_VGetArrayPointer((N_Vector)m_yp);
    m_y_out .insert(m_y_out .end(), ydata,  ydata  + neq);
    m_yp_out.insert(m_yp_out.end(), ypdata, ypdata + neq);
}

// ---------------------------------------------------------------------------
// HDF5 streaming — platform-level implementation.
// ---------------------------------------------------------------------------

void FredSolverBase::openH5File(const std::string& app_name)
{
    if (m_output_file.empty()) return;

    const int nz  = m_geom.nz;
    const int nf  = m_geom.nf;
    const int nc  = m_geom.nc;
    const int nT  = nz * (nf + nc);
    const int neq = getSolverNeq();

    hid_t file = H5Fcreate(m_output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    m_h5ctx.file     = file;
    m_h5ctx.nsteps   = 0;
    m_h5ctx.nT_cols  = nT;
    m_h5ctx.neq_cols = neq;

    // /metadata
    hid_t meta = H5Gcreate(file, "metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5_attr_str(meta, "app", app_name.c_str());
    h5_attr_int(meta, "nz",  nz);
    h5_attr_int(meta, "nf",  nf);
    h5_attr_int(meta, "nc",  nc);
    h5_attr_int(meta, "n_steps", 0);
    H5Gclose(meta);

    // /time
    m_h5ctx.ds_time = h5_create_1d(file, "time");

    // /thermal  (kept open so apps can add datasets to it)
    hid_t therm = H5Gcreate(file, "thermal", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_h5ctx.grp_therm = therm;
    m_h5ctx.ds_T      = h5_create_2d(therm, "T",           (hsize_t)nT);
    m_h5ctx.ds_peak_T = h5_create_1d(therm, "peak_T_fuel");

    // /restart
    hid_t rst = H5Gcreate(file, "restart", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_h5ctx.ds_y  = h5_create_2d(rst, "y",  (hsize_t)neq, 1);
    m_h5ctx.ds_yp = h5_create_2d(rst, "yp", (hsize_t)neq, 1);
    H5Gclose(rst);

    // App-specific datasets (mechanical/, burnup/, etc.)
    openAppH5Datasets();
}

void FredSolverBase::closeH5File()
{
    if (m_h5ctx.file < 0) return;

    // Update n_steps attribute to final count.
    hid_t meta = H5Gopen((hid_t)m_h5ctx.file, "metadata", H5P_DEFAULT);
    h5_attr_int_update(meta, "n_steps", (int)m_h5ctx.nsteps);
    H5Gclose(meta);

    closeAppH5Datasets();

    if (m_h5ctx.ds_time  >= 0) H5Dclose((hid_t)m_h5ctx.ds_time);
    if (m_h5ctx.ds_T     >= 0) H5Dclose((hid_t)m_h5ctx.ds_T);
    if (m_h5ctx.ds_peak_T>= 0) H5Dclose((hid_t)m_h5ctx.ds_peak_T);
    if (m_h5ctx.ds_y     >= 0) H5Dclose((hid_t)m_h5ctx.ds_y);
    if (m_h5ctx.ds_yp    >= 0) H5Dclose((hid_t)m_h5ctx.ds_yp);
    if (m_h5ctx.grp_therm>= 0) H5Gclose((hid_t)m_h5ctx.grp_therm);
    H5Fclose((hid_t)m_h5ctx.file);

    m_h5ctx = FredH5Ctx{};
}

void FredSolverBase::flushH5Step(double t)
{
    if (m_h5ctx.file < 0) return;

    const hsize_t step  = (hsize_t)m_h5ctx.nsteps;
    const int     nT    = m_h5ctx.nT_cols;
    const int     neq   = m_h5ctx.neq_cols;
    const int     nf    = m_geom.nf;
    const int     nc    = m_geom.nc;
    const int     nz    = m_geom.nz;

    // /time
    h5_append_scalar((hid_t)m_h5ctx.ds_time, step, t);

    // /thermal/T  — last nT entries of m_T_out (just pushed by storeOutput)
    if (nT > 0 && (int)m_T_out.size() >= nT) {
        const double* T_ptr = m_T_out.data() + m_T_out.size() - nT;
        h5_append_row((hid_t)m_h5ctx.ds_T, step, T_ptr, (hsize_t)nT);

        // /thermal/peak_T_fuel — max over all fuel nodes in the current step
        double peak = 0.0;
        for (int j = 0; j < nz; ++j)
            for (int i = 0; i < nf; ++i)
                peak = std::max(peak, T_ptr[j*(nf+nc) + i]);
        h5_append_scalar((hid_t)m_h5ctx.ds_peak_T, step, peak);
    }

    // /restart/y, /restart/yp — last neq entries of m_y_out / m_yp_out
    if (neq > 0 && (int)m_y_out.size() >= neq) {
        const double* y_ptr  = m_y_out .data() + m_y_out .size() - neq;
        const double* yp_ptr = m_yp_out.data() + m_yp_out.size() - neq;
        h5_append_row((hid_t)m_h5ctx.ds_y,  step, y_ptr,  (hsize_t)neq);
        h5_append_row((hid_t)m_h5ctx.ds_yp, step, yp_ptr, (hsize_t)neq);
    }

    // App scalars (gap_width, gas_pressure, etc.)
    appendAppH5Row();

    ++m_h5ctx.nsteps;

    // Flush HDF5 library buffers to OS after every output step so that
    // partial results are readable if the run is interrupted.
    H5Fflush((hid_t)m_h5ctx.file, H5F_SCOPE_GLOBAL);

    // Trim all in-memory output to bounded size (keep only latest step).
    if (!m_times_out.empty()) {
        double last_t = m_times_out.back();
        m_times_out = {last_t};
    }
    if (nT > 0 && (int)m_T_out.size() > nT) {
        m_T_out = std::vector<double>(m_T_out.end() - nT, m_T_out.end());
    }
    if (neq > 0 && (int)m_y_out.size() > neq) {
        m_y_out  = std::vector<double>(m_y_out .end() - neq, m_y_out .end());
        m_yp_out = std::vector<double>(m_yp_out.end() - neq, m_yp_out.end());
    }
    trimAppOutputVectors();
}

// ---------------------------------------------------------------------------
// logStepOutput — default per-step log line.
// ---------------------------------------------------------------------------
void FredSolverBase::logStepOutput(double tret, double dt_next) {
    std::cout << "  t = " << tret << " s   next dt = " << dt_next << " s\n";
}

// ---------------------------------------------------------------------------
// writeFredState — binary serialiser for checkpoint and snapshot (member).
// ---------------------------------------------------------------------------
void FredSolverBase::writeFredState(std::ofstream& f, double t) {
    const int neq = getSolverNeq();
    const int nz  = m_geom.nz;

    f.write("FREDCKPT", 8);
    uint32_t ver = 1;
    f.write(reinterpret_cast<const char*>(&ver), 4);
    f.write(reinterpret_cast<const char*>(&t),   8);

    int32_t ineq = static_cast<int32_t>(neq);
    int32_t inz  = static_cast<int32_t>(nz);
    f.write(reinterpret_cast<const char*>(&ineq), 4);
    f.write(reinterpret_cast<const char*>(&inz),  4);

    const double* y  = N_VGetArrayPointer((N_Vector)m_y);
    const double* yp = N_VGetArrayPointer((N_Vector)m_yp);
    f.write(reinterpret_cast<const char*>(y),  neq * 8);
    f.write(reinterpret_cast<const char*>(yp), neq * 8);

    RodResiduals& res = rodResiduals();
    for (int j = 0; j < nz; ++j) {
        int32_t gopen = res.isGapOpen(j) ? 1 : 0;
        f.write(reinterpret_cast<const char*>(&gopen), 4);
    }
    for (int j = 0; j < nz; ++j) {
        double off = res.axialOffset(j);
        f.write(reinterpret_cast<const char*>(&off), 8);
    }

    writeAppCheckpoint(f);
}

// ---------------------------------------------------------------------------
// readFredState — binary deserialiser (member).  Returns time from file.
// ---------------------------------------------------------------------------
double FredSolverBase::readFredState(std::ifstream& f, const std::string& filename) {
    char magic[8];
    f.read(magic, 8);
    if (std::string(magic, 8) != "FREDCKPT")
        throw std::runtime_error("loadCheckpoint: invalid magic in " + filename);

    uint32_t ver;
    f.read(reinterpret_cast<char*>(&ver), 4);
    if (ver != 1)
        throw std::runtime_error("loadCheckpoint: unsupported version");

    double t;
    f.read(reinterpret_cast<char*>(&t), 8);

    int32_t neq, nz;
    f.read(reinterpret_cast<char*>(&neq), 4);
    f.read(reinterpret_cast<char*>(&nz),  4);

    m_restart_y .resize(neq);
    m_restart_yp.resize(neq);
    f.read(reinterpret_cast<char*>(m_restart_y .data()), neq * 8);
    f.read(reinterpret_cast<char*>(m_restart_yp.data()), neq * 8);

    m_restart_gap_open    .resize(nz);
    m_restart_axial_offset.resize(nz);
    for (int j = 0; j < nz; ++j) {
        int32_t g;
        f.read(reinterpret_cast<char*>(&g), 4);
        m_restart_gap_open[j] = (g != 0);
    }
    f.read(reinterpret_cast<char*>(m_restart_axial_offset.data()), nz * 8);

    readAppCheckpoint(f);
    return t;
}

// ---------------------------------------------------------------------------
// saveCheckpoint — overwrite single fault-recovery file.
// ---------------------------------------------------------------------------
void FredSolverBase::saveCheckpoint(double t) {
    std::string fname = m_ckpt_prefix + "_checkpoint.chk";
    std::ofstream f(fname, std::ios::binary | std::ios::trunc);
    if (!f) { std::cerr << "saveCheckpoint: cannot open " << fname << "\n"; return; }

    writeFredState(f, t);
    std::cout << "  [checkpoint: " << fname << "  t=" << t << " s]\n";
}

// ---------------------------------------------------------------------------
// saveSnapshot — write a permanent restart file.
// ---------------------------------------------------------------------------
void FredSolverBase::saveSnapshot(double t, const std::string& filename) {
    std::ofstream f(filename, std::ios::binary | std::ios::trunc);
    if (!f) { std::cerr << "saveSnapshot: cannot open " << filename << "\n"; return; }

    writeFredState(f, t);
    std::cout << "  [snapshot saved: " << filename << "  t=" << t << " s]\n";
}

// ---------------------------------------------------------------------------
// setSnapshotPrefix — enable snapshot saving for the next run().
// ---------------------------------------------------------------------------
void FredSolverBase::setSnapshotPrefix(const std::string& prefix,
                                        std::vector<double> snapshot_times) {
    m_snapshot_prefix  = prefix;
    m_snapshot_timings = std::move(snapshot_times);
}

// ---------------------------------------------------------------------------
// loadCheckpoint — restore state at saved time (fault recovery).
// ---------------------------------------------------------------------------
void FredSolverBase::loadCheckpoint(const std::string& filename) {
    std::ifstream f(filename, std::ios::binary);
    if (!f) throw std::runtime_error("loadCheckpoint: cannot open " + filename);

    double t = readFredState(f, filename);
    m_restart_time        = t;
    m_restart_is_snapshot = false;
    std::cout << "loadCheckpoint: t=" << t << " s  (" << filename << ")\n";
}

// ---------------------------------------------------------------------------
// loadSnapshot — restore physical state; simulation time resets to 0.0.
// ---------------------------------------------------------------------------
void FredSolverBase::loadSnapshot(const std::string& filename) {
    std::ifstream f(filename, std::ios::binary);
    if (!f) throw std::runtime_error("loadSnapshot: cannot open " + filename);

    double t = readFredState(f, filename);
    m_restart_time        = 0.0;   // time resets for the new simulation
    m_restart_is_snapshot = true;
    std::cout << "loadSnapshot: physical state at t=" << t
              << " s; new simulation starts at t=0  (" << filename << ")\n";
}

} // namespace fred
