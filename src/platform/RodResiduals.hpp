#pragma once
#include "FuelRodGeometry.hpp"
#include "HeatConduction.hpp"
#include "StressStrain.hpp"
#include "GapStateManager.hpp"
#include "SubchannelMode.hpp"
#include "Constants.hpp"
#include <vector>

// Forward-declare SUNDIALS types.
typedef struct _generic_N_Vector *N_Vector;

namespace fred {

class FuelPelletMaterial;
class CladdingMaterial;

// \coderef{rod_residuals_class}
// -----------------------------------------------------------------------
// RodResiduals — shared infrastructure for all fuel-rod DAE assemblers.
//
// Owns the platform building blocks (HeatConduction, StressStrain,
// GapStateManager) and drives the shared residual loop:
//
//   computeResiduals()  ← NVI: calls the virtual hooks below
//       unpackState()              — fill typed state from flat IDA vector
//       computeGlobalResiduals()   — global AEs / ODEs (pressure, gas, …)
//       per-layer loop:
//           prepareLayer()             — hook for per-layer pre-work
//           computeThermalResiduals()  — thermal block for layer j
//           computeIrradiationResiduals() — irradiation block (default: no-op)
//           mechanics (in base): m_mech.computeResiduals(layerState(j), …)
//
// Subclasses supply:
//   • Typed state struct + pack/unpack
//   • Global equation count (globalOffset) and pressure accessor
//   • computeGlobalResiduals, computeThermalResiduals
//   • Optionally: computeIrradiationResiduals, freezeThermalResiduals, prepareLayer
//   • initAlgebraicState / solveMechanicalIC / setGapClosed (type-specific)
//   • gapRoot (app-specific root-finding)
// -----------------------------------------------------------------------
class RodResiduals {
public:
    // ---- Boundary conditions ----
    void setPowerDensity(TimeSeries fn)                          { m_powerFn         = std::move(fn); }
    void setLayerPowerFunctions  (std::vector<TimeSeries> fns)  { m_layerPowerFns   = std::move(fns); }
    void setCoolantTemperature(TimeSeries fn)                    { m_coolantTFn      = std::move(fn); }
    void setLayerCoolantFunctions(std::vector<TimeSeries> fns)  { m_layerCoolantFns = std::move(fns); }
    void setCoolantPressure(double p)                            { m_pcool           = p; }
    void setCoolantPressure(TimeSeries fn)                       { m_pcoolFn         = std::move(fn); }
    virtual void setInitialTemperature(double T0)                { m_T_init          = T0; }

    // ---- Physics toggles ----
    void setEnableHeatConduction(bool b) { m_heat_on = b; }
    void setEnableStressStrain(bool b)   { m_mech_on = b; }
    bool isHeatOn()              const   { return m_heat_on; }
    bool isMechOn()              const   { return m_mech_on; }

    // Irradiation physics on/off (burnup, fission gas, swelling, creep, and
    // — for apps whose per-step physics runs outside the DAE residual, e.g.
    // FredMNaSolver::afterAcceptedStep — GRSIS/Zr redistribution/cladding
    // wastage). Used by FredSolverBase::runHotStart() to freeze irradiation
    // while marching the thermal/mechanical state to steady state; default
    // on for normal transient runs.
    void setEnableIrradiation(bool b) { m_irr_on = b; }
    bool isIrradiationOn()      const { return m_irr_on; }

    // During FredSolverBase::runHotStart()'s pseudo-time march, boundary
    // conditions (power, coolant temperature/pressure, plenum temperature)
    // must stay pinned at t_fixed even though the underlying integrator's
    // own t advances through the pseudo-time horizon. Not used elsewhere.
    void setBcTimeOverride(bool on, double t_fixed = 0.0) {
        m_bc_time_override = on; m_bc_time_fixed = t_fixed;
    }

    // ---- Gap state ----
    virtual void setGapOpen(int layer, bool open) { m_gapMgr.setGapOpen(layer, open); }
    bool isGapOpen(int layer)  const      { return m_gapMgr.isGapOpen(layer); }

    // Record gap closure from the current IDA y/yp arrays and update gap manager.
    // NVI: calls unpackState() then onGapClosed().
    void setGapClosed(int layer, const double* y_data, const double* yp_data);

    // Checkpoint restore: set gap open/closed + axial offset for a layer.
    void restoreGapState(int layer, bool open, double axialOffset) {
        m_gapMgr.setGapOpen(layer, open);
        m_gapMgr.setAxialOffset(layer, axialOffset);
    }

    // Expose axial offset for checkpoint saving.
    double axialOffset(int layer) const { return m_gapMgr.axialOffset(layer); }

    // ---- Subchannel coolant mode ----
    void setSubchannelMode(SubchannelMode& mode) { m_subchannel = &mode; }

    int neqPerLayer() const { return m_neq_j; }

    // ---- Per-layer parallelism ----
    // Number of OpenMP threads to use for the per-layer loop in
    // computeResiduals() (see RodResiduals.cpp). 1 = serial (default,
    // matches pre-threading behavior exactly). Clamped to 1 if any bound
    // material is not thread-safe (see FuelPelletMaterial::isThreadSafe()) —
    // m_materials_thread_safe is set once at construction from the concrete
    // fuel/clad/gap materials actually passed in, since that's the only
    // place RodResiduals sees them by their most-derived type (Python
    // subclass trampolines included).
    void setNumThreads(int n) {
        m_num_threads = (n > 0 ? n : 1);
        if (m_num_threads > 1 && !m_materials_thread_safe) {
            m_num_threads = 1;
        }
    }
    int numThreads() const { return m_num_threads; }

    // ---- Shared IDA residual entry point (NVI) ----
    // \coderef{rod_residuals_nvi}
    int computeResiduals(double t, const double* y, const double* yp, double* r);

    // IDA callbacks — cast user_data to RodResiduals* and dispatch.
    static int idaResidual(double t, N_Vector y, N_Vector yp, N_Vector r, void* user_data);
    static int gapRoot    (double t, N_Vector y, N_Vector yp, double* gout, void* user_data);

protected:
    // neq_irr: additional per-layer equations for the irradiation block.
    //   0 for FRED-ROD (no irradiation), 1+nf for FRED-OX (burnup + swelling).
    RodResiduals(const FuelRodGeometry&    geom,
                 const FuelPelletMaterial& fuel,
                 const CladdingMaterial&   clad,
                 const GapMaterial&        gap_mat,
                 double                    coolant_pressure_MPa,
                 int                       neq_irr = 0);

    virtual ~RodResiduals() = default;

    // ---- Per-layer BC lookup helpers ----
    double layerPower  (int j, double t) const;
    double layerCoolant(int j, double t) const;
    double layerPcool  (double t)        const;
    double layerHTC    (int j)           const;

    // Remaps t to the frozen hot-start time when setBcTimeOverride(true, ...)
    // is active; identity otherwise. computeResiduals() applies this once,
    // centrally, to the t it hands to computeGlobalResiduals/
    // computeThermalResiduals; subclasses that look up their own extra
    // time-dependent BCs (e.g. plenum temperature) from the t they are
    // given automatically inherit the override with no code of their own.
    double bcTime(double t) const { return m_bc_time_override ? m_bc_time_fixed : t; }

    // ---- Virtual hooks consumed by computeResiduals ----
    // \coderef{rod_residuals_hooks}
    // Unpack the flat IDA y/yp into the application's typed state struct.
    virtual void unpackState(const double* y, const double* yp) const = 0;

    // Number of global equations at the front of the y-vector (1 for ROD, 3 for OX).
    virtual int globalOffset() const = 0;

    // Gas pressure from the unpacked state, used in the mechanics block.
    virtual double globalPressure() const = 0;

    // Per-layer state reference, passed to StressStrain.
    virtual const AxialLayerState& layerState(int j) const = 0;

    // Global residuals (pressure AE, gas ODEs, …).
    virtual void computeGlobalResiduals(double t, const double* yp, double* r) const = 0;

    // Thermal residuals for layer j.
    virtual void computeThermalResiduals(int j, double t,
                                          const double* yj, const double* ypj,
                                          double* rj) const = 0;

    // Thermal residuals when heat conduction is disabled.
    // Default: zero all neq_th entries in rj.  ROD overrides to pin AEs explicitly.
    virtual void freezeThermalResiduals(int j, const double* yj, double* rj) const;

    // Irradiation residuals for layer j (burnup, swelling, …).
    // Default: no-op (FRED-ROD has no irradiation block).
    virtual void computeIrradiationResiduals(int j,
                                              const double* ypj,
                                              double* rj) const {}

    // Irradiation residuals when irradiation is disabled (hot start).
    // Default: pin every irradiation-block ODE slot to yp=0 (freeze at
    // whatever value it currently holds) — correct as long as the whole
    // irradiation block is differential (true for FRED-OX and FRED-M-Na's
    // per-layer bup/efs/ec block; FRED-ROD's neq_irr is 0, so this is a
    // no-op there regardless). Offset by m_neq_th, same as
    // computeIrradiationResiduals: rj is the full per-layer residual block
    // (thermal first, then irradiation, then mechanics), not a block
    // private to this hook.
    virtual void freezeIrradiationResiduals(int /*j*/, const double* ypj, double* rj) const {
        for (int k = 0; k < m_neq_irr; ++k) rj[m_neq_th + k] = ypj[m_neq_th + k];
    }

    // Called at the start of each layer in the residual loop.
    // Default: no-op.  FRED-OX overrides to refresh s.gapOpen from the gap manager.
    virtual void prepareLayer(int j) const {}

    // Called once, serially, before the (possibly parallel) per-layer loop
    // in computeResiduals(), with the thread count about to be used.
    // Default: no-op. Subclasses whose per-layer thermal/mechanical hooks
    // mutate shared, non-const material objects (e.g. FRED-OX's
    // FredOxMOX::setBurnup()/FredOxGapMaterial::setLinearPower(), which are
    // read back inside the same call) override this to (re)size one clone
    // of that mutable material state per worker thread, so concurrent
    // layers never contaminate each other's burnup/linear-power inputs.
    virtual void prepareThreadLocalState(int /*nthreads*/) const {}

    // Called by setGapClosed() after unpackState() to record the closure offset
    // and synchronise the typed layer state.
    virtual void onGapClosed(int layer) = 0;

    // ---- Platform blocks ----
    const FuelRodGeometry&  m_geom;
    const CladdingMaterial& m_clad;
    HeatConduction          m_heat;
    StressStrain            m_mech;

    // ---- Per-layer equation counts ----
    // \coderef{rod_residuals_neq_layout}
    int m_neq_th;   // thermal block     (= HeatConduction::neq())
    int m_neq_irr;  // irradiation block (0 for FRED-ROD)
    int m_neq_mech; // mechanics block   (= StressStrain::neq())
    int m_neq_j;    // total per layer   (= m_neq_th + m_neq_irr + m_neq_mech)

    // ---- Boundary conditions ----
    double     m_pcool  = 0.0;
    double     m_T_init = T_REF;
    TimeSeries m_powerFn, m_coolantTFn, m_pcoolFn;
    std::vector<TimeSeries> m_layerPowerFns;    // per-layer override; empty = use m_powerFn
    std::vector<TimeSeries> m_layerCoolantFns;  // per-layer override; empty = use m_coolantTFn

    // ---- Physics toggles ----
    bool m_heat_on = true;
    bool m_mech_on = true;
    bool m_irr_on  = true;

    // ---- Hot-start BC-time override (see setBcTimeOverride) ----
    bool   m_bc_time_override = false;
    double m_bc_time_fixed    = 0.0;

    // ---- Subchannel mode (nullptr = prescribed Dirichlet BC) ----
    SubchannelMode* m_subchannel = nullptr;

    // ---- Gap state ----
    GapStateManager m_gapMgr;

    // ---- Per-layer parallelism state (setter/getter are public, above) ----
    int  m_num_threads = 1;
    bool m_materials_thread_safe = true;

    // ---- Work arrays ----
    // Thread-local scratch buffers, one per OpenMP worker thread, so
    // concurrent layers in the parallel-for loop never tear each other's
    // scratch data (each thread gets its own copy via `thread_local`; sizing
    // via resize() is a no-op once the first call per thread has run).
    // Replaces the old single mutable-member buffers (m_r_th/m_r_mech),
    // which were a data race once the per-layer loop is parallelized.
    std::vector<double>& threadLocalThermalBuffer() const {
        thread_local std::vector<double> buf;
        buf.resize(m_neq_th);
        return buf;
    }
    std::vector<double>& threadLocalMechBuffer() const {
        thread_local std::vector<double> buf;
        buf.resize(m_neq_mech);
        return buf;
    }
};

} // namespace fred
