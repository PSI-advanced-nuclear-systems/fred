#pragma once
#include "FuelRodGeometry.hpp"
#include "TimeTable.hpp"
#include "Constants.hpp"
#include <cstdint>
#include <vector>
#include <string>
#include <iosfwd>

// Forward-declare SUNDIALS types so IDA headers stay out of this header.
typedef struct _generic_N_Vector *N_Vector;
typedef int (*IDAResFn) (double, N_Vector, N_Vector, N_Vector, void*);
typedef int (*IDARootFn)(double, N_Vector, N_Vector, double*, void*);

namespace fred {

class RodResiduals;

// HDF5 streaming state — all handles stored as int64_t (hid_t = int64_t in
// HDF5 1.12+) to keep hdf5.h out of this header.
struct FredH5Ctx {
    int64_t file       = -1;
    int64_t grp_therm  = -1;  // thermal/ group kept open; apps add datasets here
    int64_t ds_time    = -1;
    int64_t ds_T       = -1;
    int64_t ds_peak_T  = -1;
    int64_t ds_y       = -1;
    int64_t ds_yp      = -1;
    int64_t nsteps     = 0;
    int     nT_cols    = 0;
    int     neq_cols   = 0;
};

// \coderef{fred_solver_base}
// -----------------------------------------------------------------------
// FredSolverBase — shared infrastructure for all FRED solver drivers,
// SUNDIALS-based or not.
//
// Eliminates duplication of:
//   • TimeTable piecewise-linear interpolation helper
//   • N_Vector/SUNContext storage and destructor cleanup (used either as
//     real SUNDIALS IDA state vectors, or — by FredMNaSolver, which does
//     not use IDA at all — as plain data containers so the checkpoint/
//     snapshot/HDF5-restart plumbing below can stay app-agnostic)
//   • Common members: geometry ref, initial temperature
//   • Common output storage: times, temperatures, y/yp state vectors
//   • Common result accessors: timePoints(), temperatures(), yOut(), ypOut()
//   • peakFuelTemperature() computation
//   • y/yp output append helper (storeYYP)
//   • Checkpoint/snapshot binary I/O and HDF5 streaming (app-agnostic)
//
// Subclasses implement one pure virtual hook used by the above:
//   getSolverNeq()       — total DAE equations
//
// Applications that drive themselves through SUNDIALS IDA (FRED-ROD,
// FRED-OX) derive from FredIdaSolverBase (see FredIdaSolverBase.hpp)
// instead of this class directly — that subclass adds the SUNDIALS-specific
// setup hooks, tolerances, and the shared runTimeLoop(). FRED-M-Na derives
// from FredSolverBase directly since it implements its own fixed-step
// backward-Euler integrator and never calls setupSundials()/runTimeLoop().
//
// Subclasses keep their own run(), initState(), storeOutput() since those
// contain app-specific logic, but call the protected helpers for the
// repetitive boilerplate.
// -----------------------------------------------------------------------
class FredSolverBase {
public:
    virtual ~FredSolverBase();

    // Verbosity level for console diagnostics during a run.
    // 0 = silent. Default 3 = base time-loop status line (current time +
    // next IDA step size) at every output step; see logStepOutput().
    // Applications may branch on higher levels in their own
    // logStepOutput()/afterAcceptedStep() overrides to print additional
    // application-specific variables of interest.
    void setVerbosity(int level) { m_verbosity = level; }
    int  verbosity() const { return m_verbosity; }

    // Stream output to an HDF5 file as the simulation runs (crash-safe).
    // Must be called before run().  Each dtout step is written immediately;
    // in-memory output vectors are trimmed to the last step only.
    void setOutputFile(const std::string& filename) { m_output_file = filename; }

    // Checkpoint API (fault recovery — single overwriting file).
    void setCheckpointPrefix(const std::string& prefix) { m_ckpt_prefix = prefix; }
    void loadCheckpoint(const std::string& filename);
    double restartTime() const { return m_restart_time; }

    // Snapshot API (permanent restart states — time resets to 0 on load).
    // snapshot_times: additional times to save snapshots; final state is always saved.
    // Files: <prefix>_frame1.snapshot, <prefix>_frame2.snapshot, ...
    void setSnapshotPrefix(const std::string& prefix,
                           std::vector<double> snapshot_times = {});
    void loadSnapshot(const std::string& filename);

    // Shared boundary-condition setters (NVI — implemented once here;
    // app-specific wiring goes through the protected hooks below).
    // Uniform power (broadcast to all layers).
    void setPowerDensityHistory   (std::vector<double> times, std::vector<double> qv);
    // Per-layer power density profiles (primary entry point).
    void setPowerDensityHistoryPerLayer(std::vector<double> times,
                                        std::vector<std::vector<double>> qqv_per_layer);
    // Uniform coolant temperature (broadcast to all layers).
    void setCoolantTemperature    (std::vector<double> times, std::vector<double> T_K);
    // Per-layer coolant temperature profiles (primary entry point).
    void setCoolantTemperaturePerLayer(std::vector<double> times,
                                       std::vector<std::vector<double>> T_per_layer);
    void setCoolantPressure       (double pcool_MPa);
    void setCoolantPressureHistory(std::vector<double> times, std::vector<double> pcool_MPa);
    void setPlenumTemperatureHistory(std::vector<double> times, std::vector<double> T_K);
    void setInitialTemperature    (double T0_K);   // sets m_T0, then calls hook

    // Subchannel coolant BC (no-op default; FredMNaSolver overrides with real implementation).
    virtual void setCoolantChannel(double /*dhyd*/, double /*xarea*/, double /*flowr*/,
                                   std::vector<double> /*T_inlet_times*/,
                                   std::vector<double> /*T_inlet_vals*/) {}

    // Result accessors common to all apps.
    const std::vector<double>& timePoints()   const { return m_times_out; }
    const std::vector<double>& temperatures() const { return m_T_out; }
    const std::vector<double>& yOut()         const { return m_y_out; }
    const std::vector<double>& ypOut()        const { return m_yp_out; }
    std::vector<double> peakFuelTemperature() const;

    virtual int neqTotal() const = 0;

    void setEnableHeatConduction(bool b);
    void setEnableStressStrain  (bool b);

    // Hot (steady-state) start. Default: off (cold start at
    // setInitialTemperature()'s T0). When enabled, run() first marches the
    // rod forward in pseudo-time at its t=0 boundary conditions (power,
    // coolant) — with irradiation physics (burnup, fission gas, swelling,
    // creep, GRSIS/Zr redistribution/cladding wastage as applicable) turned
    // off — until the temperature/stress-strain state stops changing, then
    // uses that converged state as t=0 for the real run instead of the cold
    // T0 state. A pre-real-time-step calculation only; must be called
    // before run().
    void setHotStart(bool b) { m_hot_start = b; }
    bool hotStart() const    { return m_hot_start; }

    // Number of OpenMP threads used for the per-axial-layer residual
    // computation during run(). 1 (default) = serial, identical behavior to
    // before threading was added. Concrete solvers forward this to their
    // RodResiduals-derived assembler (setNumThreads()) at the start of
    // run(); some materials (e.g. Python-subclassed FRED-ROD materials) are
    // not thread-safe and may clamp this back down to 1 internally.
    void setThreads(int n) { m_threads = n > 0 ? n : 1; }
    int  threads() const   { return m_threads; }

protected:
    explicit FredSolverBase(const FuelRodGeometry& geom);

    // ---- Physics toggle accessor ----
    virtual RodResiduals& rodResiduals() = 0;

    // Runs the hot-start pseudo-time march described above. Pure virtual so
    // every concrete solver is forced to define what "steady state" means
    // for its own time-integration scheme (FredIdaSolverBase provides one
    // shared IDA-based implementation for FRED-ROD/FRED-OX; FRED-M-Na's
    // fixed-step backward-Euler integrator implements its own). Only called
    // by run() when m_hot_start is true and the run is not a checkpoint/
    // snapshot restart.
    virtual void runHotStart() = 0;

    // ---- Boundary-condition hooks (called by NVI setters above) ----
    virtual void onPowerDensityPerLayerSet  (const std::vector<TimeTable>& tabs) = 0;
    virtual void onCoolantTemperaturePerLayerSet(const std::vector<TimeTable>& tabs) = 0;
    virtual void onCoolantPressureSet   (const TimeTable& tbl) = 0;
    virtual void onInitialTemperatureSet(double T0)            = 0;
    // Optional: default no-op (only meaningful for gas-bond models).
    virtual void onPlenumTemperatureSet (const TimeTable& /*tbl*/) {}

    // ---- Subclass hook shared by checkpoint/HDF5 plumbing ----
    virtual int    getSolverNeq()                           const = 0;

    // Append current IDA y/yp pointers to m_y_out / m_yp_out.
    void storeYYP();

    // Checkpoint / snapshot helpers.
    void saveCheckpoint(double t);                                    // overwrites single .chk file
    void saveSnapshot(double t, const std::string& filename);        // writes permanent .snapshot file

    // Binary serialisation helpers (private implementation detail).
    void   writeFredState(std::ofstream& f, double t);
    double readFredState (std::ifstream& f, const std::string& filename);

    // ---- Virtual hooks shared by both time-integration schemes ----

    // Unpack current y/yp into the app's state struct.
    virtual void unpackCurrentState() = 0;

    // Append a time-point snapshot to the app's output vectors.
    virtual void storeOutput(double t) = 0;

    // Called once per accepted step (IDA's runTimeLoop, or FredMNaSolver's
    // own runOneStepLoop). Default: no-op.  FredMNaSolver overrides to run
    // per-step physics (GRSIS bubble evolution, burnup, Zr redistribution,
    // cladding wastage).
    virtual void afterAcceptedStep(double /*t*/, double /*dt*/) {}

    // Per-output-step log line, called only when m_verbosity >= 3.
    // Default: current time + next step size (dt_next).
    // OX/M-Na override to also print gpres / fgrel / etc.
    virtual void logStepOutput(double tret, double dt_next);

    // Checkpoint virtual hooks — default no-ops (base handles y/yp/gap state).
    virtual void writeAppCheckpoint(std::ostream& /*os*/) const {}
    virtual void readAppCheckpoint (std::istream& /*is*/) {}

    // ---- HDF5 streaming helpers ----
    // openH5File: create file + base datasets; call at start of run() before first storeOutput.
    void openH5File(const std::string& app_name);
    // closeH5File: flush, update n_steps metadata, close all handles.
    void closeH5File();
    // flushH5Step: write current step to file then trim in-memory vectors.
    // Call at the end of each app's storeOutput() (after storeYYP()).
    void flushH5Step(double t);

    // App-specific H5 hooks (default: no-op).
    virtual void openAppH5Datasets   () {}  // create /mechanical, /burnup datasets etc.
    virtual void appendAppH5Row      () {}  // write last entry of each app scalar vector
    virtual void trimAppOutputVectors() {}  // keep only last entry of each app scalar vector
    virtual void closeAppH5Datasets  () {}  // close app dataset/group handles

    // ---- Shared members ----
    const FuelRodGeometry& m_geom;
    std::string  m_output_file;
    FredH5Ctx    m_h5ctx;
    double m_T0   = T_REF;
    int    m_verbosity = 3;
    bool   m_hot_start = false;
    int    m_threads   = 1;

    // N_Vector/SUNContext (opaque pointers keep SUNDIALS headers out of this
    // header). Used as real IDA state vectors by FredIdaSolverBase-derived
    // apps, or as plain data containers by FredMNaSolver (see class comment).
    void* m_sunctx  = nullptr;
    void* m_y       = nullptr;
    void* m_yp      = nullptr;

    // Common output storage.
    std::vector<double> m_times_out;
    std::vector<double> m_T_out;
    std::vector<double> m_y_out;
    std::vector<double> m_yp_out;

    // Checkpoint state (fault recovery).
    std::string         m_ckpt_prefix  = "";
    double              m_restart_time = -1.0;   // <0 = no restart
    bool                m_restart_is_snapshot = false;  // true when loaded via loadSnapshot
    std::vector<double> m_restart_y, m_restart_yp;
    std::vector<bool>   m_restart_gap_open;
    std::vector<double> m_restart_axial_offset;

    // Snapshot state (permanent restart states).
    std::string         m_snapshot_prefix   = "";
    std::vector<double> m_snapshot_timings  = {};
    int                 m_snapshot_count    = 0;
};

} // namespace fred
