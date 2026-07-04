#pragma once
#include "FredSolverBase.hpp"

namespace fred {

// \coderef{fred_ida_solver_base}
// -----------------------------------------------------------------------
// FredIdaSolverBase — SUNDIALS IDA-specific solver infrastructure.
//
// Split out of FredSolverBase so that applications which do not use
// SUNDIALS IDA (FRED-M-Na implements its own fixed-step backward-Euler
// "one-step" integrator — see FredMNaSolver's class comment) are not forced
// to implement an IDA residual/Jacobian entry point they never call.
// FRED-ROD and FRED-OX — the two apps that do drive themselves through
// IDASolve — derive from this class instead of FredSolverBase directly.
//
// Adds on top of FredSolverBase:
//   • rtol/atol/hmax/hinit/max_nonlin_iters/max_ord SUNDIALS tuning knobs
//   • SUNDIALS IDA object storage (m_ida_mem, m_A, m_LS) and its destructor
//     cleanup (m_y/m_yp/m_sunctx remain in FredSolverBase — they are also
//     used by FredMNaSolver as plain N_Vector data containers)
//   • setupSundials() / setupGapRoots() / calcIC() / reinitAfterEvent() /
//     applyRestartToIDA() SUNDIALS helpers
//   • runTimeLoop() — the shared IDA one-step time-integration loop, and
//     the gap-event dispatch hooks (isGapOpen/handleGapClosed/
//     handleGapReopened/afterGapEventsHandled) it alone calls
//
// Subclasses implement four pure virtual hooks consumed by setupSundials():
//   getSolverResFn()     — IDA residual callback
//   getSolverUserData()  — pointer passed to IDA (usually &m_res)
//   fillSolverIdArray()  — fill differential/algebraic id vector
//   packSolverState()    — pack application state into the IDA y vector
// (getSolverNeq() — the fifth hook from the original undivided class — stays
// pure virtual on FredSolverBase itself: it is also needed by the
// checkpoint/HDF5 plumbing shared with FredMNaSolver.)
// -----------------------------------------------------------------------
class FredIdaSolverBase : public FredSolverBase {
public:
    ~FredIdaSolverBase() override;

    void setTolerances     (double rtol, double atol) { m_rtol = rtol; m_atol = atol; }
    void setMaxStep        (double hmax)    { m_hmax = hmax; }
    void setInitStep       (double hinit)   { m_hinit = hinit; }
    void setMaxNonlinIters (int n)          { m_max_nonlin_iters = n; }
    // Maximum BDF order (IDASetMaxOrd). IDA's BDF default max order is 5.
    // Forcing a lower order (e.g. 1 = backward Euler) trades local
    // truncation-error accuracy for a smoother, less stiff nonlinear solve.
    // Valid range 1-5; <=0 means "not set" (SUNDIALS default applies).
    void setMaxOrd         (int order)      { m_max_ord = order; }

protected:
    explicit FredIdaSolverBase(const FuelRodGeometry& geom);

    // ---- Subclass hooks consumed by setupSundials() ----
    virtual IDAResFn getSolverResFn()                      const = 0;
    virtual void*  getSolverUserData()                           = 0;
    virtual void   fillSolverIdArray(double* id, int neq)  const = 0;
    virtual void   packSolverState  (double* y)            const = 0;

    // ---- SUNDIALS helpers ----

    // Create SUNDIALS context, allocate N_Vectors, initialise IDA, attach
    // dense linear solver, set tolerances, mark id array, suppress algebraic
    // variables from error test.  Calls packSolverState to fill y.
    void setupSundials(int max_steps = 10000, double t0 = 0.0);

    // Register gap root detection.  Call after setupSundials().
    void setupGapRoots(int nz, IDARootFn rootFn);

    // IDACalcIC with consistent-IC fallback.  Returns true on success.
    bool calcIC(double t_ic);

    // IDAReInit + IDACalcIC with save/restore on CalcIC failure.
    // Call this inside the time loop after all gap-state changes are applied.
    void reinitAfterEvent(double tret, double dtout);

    // Shared IDA time-integration loop.  Call after setupSundials(),
    // optional setupGapRoots(), calcIC(), and an initial storeOutput(0.0).
    // Dispatches into the virtual hooks below on each step and gap event.
    void runTimeLoop(double tend, double dtout, bool all_steps, double t_start = 0.0);

    // Inject checkpoint/restart state into a freshly initialised IDA.
    void applyRestartToIDA(double t_start);

    // FredSolverBase::runHotStart() — shared IDA-based implementation for
    // FRED-ROD/FRED-OX. Reuses the IDA object already built by the caller's
    // setupSundials()+calcIC() (must be called first): marches forward in
    // pseudo-time at the current (t=0) boundary conditions, with
    // irradiation physics disabled and boundary-condition lookups pinned to
    // t=0 (RodResiduals::setBcTimeOverride), via repeated IDASolve(...,
    // IDA_NORMAL) calls to a geometrically-growing pseudo-time target, until
    // the differential (id==1) part of y stops changing. Re-anchors IDA at
    // t=0 with the converged state via IDAReInit + calcIC() before
    // returning, so the caller's subsequent real time loop starts fresh.
    void runHotStart() override;

    // ---- Gap-state queries / mutations, used only by runTimeLoop's ----
    // ---- IDA root-event dispatch.                                  ----
    virtual bool isGapOpen     (int j);
    virtual void handleGapClosed  (int j);
    virtual void handleGapReopened(int j);

    // Called once after all gap events at tret are applied, before returning
    // to the IDA loop.  Default: call reinitAfterEvent(tret, dtout).
    // ROD overrides to also re-solve mechanical equilibrium.
    virtual void afterGapEventsHandled(double tret, double dtout);

    // ---- SUNDIALS tuning knobs ----
    double m_rtol = 1.0e-6;
    double m_atol = 1.0e-8;
    double m_hmax = -1.0;          // IDASetMaxStep   (<= 0 → not set)
    double m_hinit = -1.0;         // IDASetInitStep  (<= 0 → not set)
    int    m_max_nonlin_iters = -1; // IDASetMaxNonlinIters (<= 0 → not set)
    int    m_max_ord = -1;         // IDASetMaxOrd    (<= 0 → not set)

    // SUNDIALS IDA objects (opaque pointers keep IDA headers out of this
    // header). m_y/m_yp/m_sunctx stay on FredSolverBase.
    void* m_ida_mem = nullptr;
    void* m_A       = nullptr;
    void* m_LS      = nullptr;
};

} // namespace fred
