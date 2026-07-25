#pragma once
#include "FredMNaResiduals.hpp"
#include "FredMNaState.hpp"
#include "FredMNaFailure.hpp"
#include "FredMNaGrsis.hpp"
#include "FredMNaCladWastage.hpp"
#include "FredMNaCladSwelling.hpp"
#include "FredMNaSodiumProperties.hpp"
#include "FredMNaSubchannelMode.hpp"
#include "fuelpelletmaterial/UPuZr.hpp"
#include "claddingmaterial/HT9.hpp"
#include "gapmaterial/Sodium.hpp"
#include "platform/FredSolverBase.hpp"
#include <memory>
#include <vector>
#include <iosfwd>

namespace fred {

// -----------------------------------------------------------------------
// FredMNaSolver — fixed-step, backward-Euler, "always accept" time
// integration driver for FRED-M-Na (its own "one-step" scheme).
//
// Unlike FRED-ROD/FRED-OX, FredMNaSolver does NOT use SUNDIALS IDA to
// advance in time. Legacy FRED-M (Baseir.for) takes a fixed time step
// dt=dtout and iterates a Picard (successive-substitution) outer loop,
// capped at a fixed iteration count, ALWAYS accepting whatever iterate it
// ends on (no step-size shrinking / retry on non-convergence, unlike IDA's
// adaptive BDF corrector). Porting that behaviour onto IDA would require
// overriding IDA's own step-acceptance logic, which IDA does not expose;
// instead FredMNaSolver implements its own fixed-step outer Picard loop
// (globals solved by direct substitution + one dense Newton solve per
// axial layer, since layers only couple through the 3 rod-scalar globals —
// see FredMNaResiduals::computeLayerResiduals/computeGlobalUpdate) that
// replicates the same "fixed dt, always accept" philosophy. See
// runOneStepLoop / solveStepBackwardEuler / newtonSolveLayer below.
//
// SUNDIALS N_Vector (m_y/m_yp, inherited from FredSolverBase) is still used
// as a plain data container so checkpoint/snapshot/HDF5-restart plumbing
// in FredSolverBase keeps working unmodified — but IDACreate/IDASolve/
// IDAReInit (the actual IDA integrator) are never called for this app.
//
// Inherits shared infrastructure from FredSolverBase and adds:
//   • Sodium (liquid bond) gap material
//   • Per-step physics: Zr redistribution, cladding wastage, GRSIS bubbles
//   • Per-layer bup_FIMA tracking outside the y-vector
//   • Failure criteria: burst margin (HT9), melt margin (UPuZr solidus)
// -----------------------------------------------------------------------
class FredMNaSolver : public FredSolverBase {
public:
    FredMNaSolver(const FuelRodGeometry&    geom,
                  UPuZr&                    fuel,
                  FredMNaCladdingMaterial&  clad);

    // Initial fill-gas pressure setter; setBondPressure is a legacy alias.
    void setInitialGasPressure(double p_MPa);
    void setBondPressure(double p_MPa) { setInitialGasPressure(p_MPa); }

    // Select thermal conductivity correction model (default: DetailedNaSodium).
    void setConductivityModel(ConductivityModel m) { m_fuel.setConductivityModel(m); }

    void setEnableZrRedistribution(bool b)  { m_enable_zr    = b; }
    void setEnableCladWastage(bool b)       { m_enable_waste = b; }
    // Select the cladding wastage model (default: PrecipitationKinetics,
    // the Clanth.for port).  LaTracking solves per-node lanthanide diffusion
    // (MFUEL App. 9.3) and integrates the interface flux.
    void setCladWastageModel(CladWastageModel m) { m_wastage_model = m; }
    // LaTracking only: gate the interface sink on soft/clos contact (default
    // true, same as the PK model) or keep it active from t=0 (false, MFUEL's
    // unconditional Dirichlet — transport through the sodium bond).
    void setLaSinkRequiresContact(bool b)
        { m_wastage_params.la_sink_requires_contact = b; }
    // LaTracking only: radial subgrid refinement factor for the La diffusion
    // solve (default 8; 1 = solve on the thermal grid, under-resolves the
    // floored-diffusivity boundary layer).
    void setLaSubgridRefine(int m)
        { m_wastage_params.la_refine = std::max(1, m); }

    // Cladding void swelling (FredMNaCladSwelling.hpp). Off by default:
    // this is new physics with no legacy-parity default run depending on
    // it, and enabling it changes gap-closure timing for every existing
    // script, so it is strictly opt-in.
    void setEnableCladSwelling(bool b)      { m_enable_swelling = b; }
    // Select the cladding void-swelling model (default: Hofmann, Cswel.for
    // icswel=1). SAS4A is Cswel.for icswel=2, gated on dose >= threshold.
    void setCladSwellingModel(CladSwellingModel m) { m_swelling_model = m; }
    // Fast-flux-to-linear-power ratio (fltpow) and fluence->dose conversion
    // (dconv) — case-specific neutronics calibration constants (Serpent-
    // derived), NOT universal physical constants. Defaults are Timpano's
    // Inner-Fuel-Average values; override for other pin/core positions.
    void setFastFluxCalibration(double fltpow, double dconv) {
        m_swelling_params.fltpow = fltpow;
        m_swelling_params.dconv  = dconv;
    }
    // SAS4A only: dose threshold [dpa] below which swelling is zero (default 100).
    void setCladSwellingDoseThreshold(double dpa)
        { m_swelling_params.dose_threshold_dpa = dpa; }
    void setEnableGrsis(bool b)             { m_enable_grsis = b; }
    void setGrsisDataMode(GrsisDataMode m)  { m_grsis_mode   = m; }
    void setSodiumMode(SodiumMode m)        { m_gap_na.~Sodium(); new (&m_gap_na) Sodium(m); }

    // Fixed internal time step for the one-step backward-Euler integrator
    // (distinct from dtout, the output-writing cadence passed to run()).
    // <=0 (default): use dtout itself as the step size (one fixed step per
    // output interval — matches legacy's dt==dtout convention when dtout is
    // set to the same cadence as the reference input deck; coarser dtout
    // values are a deliberate performance/accuracy trade-off left to the
    // caller, since backward Euler is unconditionally stable — see class
    // comment). A step is always shortened to avoid overshooting the next
    // output time, so a larger step_size than dtout has no effect.
    void   setStepSize(double dt) { m_step_size = dt; }
    double stepSize()       const { return m_step_size; }

    // Subchannel coolant BC: computes T_co axially from inlet + mass-energy balance;
    // evaluates cladding-to-coolant HTC from Peclet-number correlation (Robin BC).
    // This is the only supported coolant BC for FredMNaSolver.
    void setCoolantChannel(double dhyd, double xarea, double flowr,
                           std::vector<double> T_inlet_times,
                           std::vector<double> T_inlet_vals,
                           FredMNaSubchannelMode::HtcCorrelation corr
                               = FredMNaSubchannelMode::HtcCorrelation::kMikityuk);

    // threads: number of OpenMP threads for the per-axial-layer Newton
    // solve in solveStepBackwardEuler (default 1 = serial). Layers only
    // couple through the 3 rod-scalar globals, which are held fixed
    // (Jacobi) for the duration of one outer sweep's per-layer loop, so
    // parallelizing across layers changes nothing about what is being
    // solved — see solveStepBackwardEuler's header comment.
    void run(double tend, double dtout, bool all_steps = false, int threads = 1);

    int neqTotal() const override { return m_res.neqTotal(); }

    // Result accessors (temperature/times/y/yp from base).
    const std::vector<double>& gasPressure() const { return m_gpres_out; }
    const std::vector<double>& fgGen()       const { return m_fggen_out; }
    const std::vector<double>& fgRel()       const { return m_fgrel_out; }
    const std::vector<double>& gapWidth()    const { return m_gap_out; }
    const std::vector<double>& burnup()      const { return m_bup_out; }
    const std::vector<double>& burnupAtPct() const { return m_bupave_atpct_out; }
    const std::vector<double>& cladWastage() const { return m_xwast_out; }

    // GRSIS swelling outputs (per time step, spatially averaged over layers)
    const std::vector<double>& grisSwellingTotal() const { return m_swtot_out; }
    const std::vector<double>& grisSwellingOpen()  const { return m_swopen_out; }

    // Failure margin outputs (per time step, worst layer)
    const std::vector<double>& burstMargin() const { return m_burst_out; }
    const std::vector<double>& meltMargin()  const { return m_melt_out; }

    // Per-layer outputs (flattened [n_times * nz]; stride = nz per time step)
    const std::vector<double>& coolantTemperaturePerLayer() const { return m_T_co_out; }
    const std::vector<double>& gapWidthPerLayer()           const { return m_gap_layer_out; }
    const std::vector<double>& hgapPerLayer()               const { return m_hgap_layer_out; }

    // Compatibility aliases for Python bindings (base class renamed these).
    const std::vector<double>& times()        const { return m_times_out; }
    const std::vector<double>& temperatures() const { return m_T_out; }
    const std::vector<double>& y_out()        const { return m_y_out; }
    const std::vector<double>& yp_out()       const { return m_yp_out; }

protected:
    // FredSolverBase boundary-condition hooks.
    RodResiduals& rodResiduals()                                             override { return m_res; }
    void onPowerDensityPerLayerSet      (const std::vector<TimeTable>& tabs) override;
    void onCoolantTemperaturePerLayerSet(const std::vector<TimeTable>& tabs) override;

    // FredMNaSolver always uses subchannel mode; prescribing coolant temperature
    // directly is not supported and will throw std::logic_error.
    void onCoolantPressureSet   (const TimeTable& tbl) override;
    void onInitialTemperatureSet(double T0)            override;
    void onPlenumTemperatureSet (const TimeTable& tbl) override;

    // FredSolverBase setup hook (still needed: shared by checkpoint/HDF5
    // plumbing). FredMNaSolver does not derive from FredIdaSolverBase (it
    // never calls setupSundials()/runTimeLoop() — see class comment above),
    // so it has no getSolverResFn()/getSolverUserData()/fillSolverIdArray()/
    // packSolverState() to implement; m_res.fillIdArray()/m_res.pack() are
    // called directly from solveStepBackwardEuler()/run() instead.
    int      getSolverNeq()                          const override { return m_res.neqTotal(); }

    // FredSolverBase runTimeLoop hooks.
    void afterAcceptedStep (double t, double dt) override;
    void unpackCurrentState()         override;
    void storeOutput       (double t) override;
    void logStepOutput     (double tret, double dt_next) override;

private:
    UPuZr&  m_fuel;
    FredMNaCladdingMaterial& m_clad;
    Sodium  m_gap_na;

    FredMNaSodiumProperties                 m_na_props;
    std::unique_ptr<FredMNaSubchannelMode>  m_subchannel_mode;

    FredMNaResiduals m_res;
    FredMNaRodState  m_state;

    std::vector<TimeTable> m_powerTab;
    std::vector<TimeTable> m_coolantTab;
    TimeTable              m_plenumTab;

    bool         m_enable_zr    = true;
    bool         m_enable_waste = true;
    CladWastageModel m_wastage_model = CladWastageModel::PrecipitationKinetics;
    CladWastageParams m_wastage_params{};
    bool m_enable_swelling = false;
    CladSwellingModel m_swelling_model = CladSwellingModel::Hofmann;
    CladSwellingParams m_swelling_params{};
    bool         m_enable_grsis = true;
    GrsisDataMode m_grsis_mode  = GrsisDataMode::FEAST;

    std::vector<double> m_gpres_out, m_fggen_out, m_fgrel_out;
    std::vector<double> m_gap_out, m_bup_out, m_bupave_atpct_out, m_xwast_out;
    std::vector<double> m_swtot_out, m_swopen_out;
    std::vector<double> m_burst_out, m_melt_out;

    // Per-layer output vectors (n_times * nz entries each)
    std::vector<double> m_T_co_out;
    std::vector<double> m_gap_layer_out;
    std::vector<double> m_hgap_layer_out;

    bool   m_grsis_first   = true;  // reset each run(); signals GRSIS init call
    double m_step_size     = -1.0;  // <=0 => use dtout as the fixed step

    // FredSolverBase::runHotStart() — FRED-M-Na's own implementation (this
    // app does not derive from FredIdaSolverBase, so it cannot reuse that
    // class's IDA-based march). Repeatedly calls solveStepBackwardEuler()
    // with a fixed BC-evaluation time of t=0 and a geometrically-growing
    // step size, with irradiation physics disabled (RodResiduals::
    // setEnableIrradiation(false)) and afterAcceptedStep() (GRSIS/Zr
    // redistribution/cladding wastage) simply not called, until the
    // differential part of y stops changing. See run().
    void runHotStart() override;

    // ---- One-step (backward-Euler, fixed dt, "always accept") integrator ----
    // Plays the role FredIdaSolverBase::runTimeLoop/setupSundials play for
    // FRED-ROD/FRED-OX (see class comment) — this app derives from
    // FredSolverBase directly, not FredIdaSolverBase. No SUNDIALS IDA
    // integrator calls anywhere below.
    void runOneStepLoop(double tend, double dtout, bool all_steps, double t_start);
    // Advances y (m_y N_Vector data) from t_new-dt to t_new in place via the
    // fixed-step outer Picard sweep (globals direct-substitution + per-layer
    // dense Newton). Always returns after MAX_OUTER sweeps or convergence —
    // never rejects/retries with a smaller dt.
    void solveStepBackwardEuler(double t_new, double dt);
    // Per-layer local Newton solve (FredMNaDenseSolve.hpp), backward-Euler
    // discretisation of layer j's own neq_j-sized block only. gpres/fggen/
    // fgrel (the 3 rod-scalar globals) are NOT unknowns here: they are held
    // fixed for the whole outer sweep (set once, serially, by
    // solveStepBackwardEuler via m_res.setGlobalState() before the
    // per-layer loop starts), so this call only ever reads them, never
    // writes them — safe to run concurrently across layers/threads. yj/ypj
    // updated in place. Returns the infinity norm of the final Newton
    // correction over this layer's own unknowns.
    double newtonSolveLayer(int j, double t_new, double dt,
                             const double* yj_old, const double* idj,
                             double* yj, double* ypj, int n) const;
    // Restart-state injection without touching IDA
    // (FredIdaSolverBase::applyRestartToIDA calls IDAReInit, which requires
    // m_ida_mem — this app has no such member, it does not derive from
    // FredIdaSolverBase).
    void applyRestartOneStep(double t_start);

    // HDF5 dataset handles for MNA-specific quantities (int64_t = hid_t).
    struct MnaH5Ctx {
        int64_t grp_burnup = -1;
        int64_t ds_gap     = -1;
        int64_t ds_gpres   = -1;
        int64_t ds_fggen   = -1;
        int64_t ds_fgrel   = -1;
        int64_t ds_bup        = -1;
        int64_t ds_bup_atpct  = -1;
        int64_t ds_xwast   = -1;
        int64_t ds_swtot   = -1;
        int64_t ds_swopen  = -1;
        int64_t ds_burst   = -1;
        int64_t ds_melt    = -1;
        // 2D (nsteps, nz*nf) per-layer/per-node datasets, flattened [j*nf+i]
        int64_t ds_kfuel     = -1;  // thermal/k_fuel   fuel k [W/(m·K)]
        int64_t ds_psod      = -1;  // thermal/psod     Na infiltration fraction [-]
        int64_t ds_ptot      = -1;  // thermal/poros_tot  porosity seen by k model [-]
        int64_t ds_pgas      = -1;  // thermal/poros_gas  gas-filled share [-]
        int64_t ds_zrwf      = -1;  // burnup/zr_wf     local Zr weight fraction [-]
        int64_t grp_grsis    = -1;  // /grsis group
        int64_t ds_sw_close  = -1;  // grsis/swclose    closed-bubble swelling [-]
        int64_t ds_sw_open   = -1;  // grsis/swopen     open-bubble swelling [-]
        int64_t ds_sw_sol    = -1;  // grsis/swsol      solid-FP swelling [-]
        int64_t ds_sw_tot    = -1;  // grsis/swtot      total swelling [-]
        int64_t ds_c_la      = -1;  // burnup/c_la      La number density [#/m3]
        // 2D (nsteps, nz) per-layer datasets
        int64_t ds_xwast_layer = -1; // burnup/xwast_layer  wastage per layer [m]
        int64_t ds_contact     = -1; // thermal/contact_state 0=open 1=soft 2=clos
        // Cladding void swelling (FredMNaCladSwelling.hpp)
        int64_t grp_swelling = -1;  // /swelling group
        int64_t ds_ecs       = -1;  // swelling/ecs (nsteps, nz*nc) linear swelling strain [-]
        int64_t ds_neuflue2  = -1;  // swelling/neuflue2 (nsteps, nz) fast fluence [1e22 n/cm2]
        int64_t ds_clad_dose = -1;  // swelling/dose (nsteps, nz) cladding dose [dpa]
    };
    MnaH5Ctx m_h5mna;

    void openAppH5Datasets   () override;
    void appendAppH5Row      () override;
    void trimAppOutputVectors() override;
    void closeAppH5Datasets  () override;

    void initState();

    // Checkpoint hooks (override base no-ops).
    void writeAppCheckpoint(std::ostream& os) const override;
    void readAppCheckpoint (std::istream& is)       override;
};

} // namespace fred
