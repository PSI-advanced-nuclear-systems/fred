#include "FredIdaSolverBase.hpp"
#include "RodResiduals.hpp"

#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_openmp.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_context.h>

#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <vector>

namespace fred {

namespace {
// SUNDIALS ships a genuinely multithreaded N_Vector (NVECTOR_OPENMP): its
// WRMS-norm, dot-product, and linear-combination kernels -- used throughout
// IDA's nonlinear solve and local error control on every step -- run across
// OpenMP threads instead of a single one. This is additive to (not a
// replacement for) the per-layer residual-assembly threading added
// elsewhere: that parallelizes computing F(y), this parallelizes the O(neq)
// vector arithmetic IDA itself performs on y/yp/the correction each step.
// The O(neq^3) dense Newton-matrix factorization (SUNLinSol_Dense) still
// runs single-threaded regardless -- SUNDIALS has no in-built threaded dense
// direct solver -- so this does not remove that bottleneck, only shrinks the
// serial vector-kernel share of per-step cost.
N_Vector makeVector(sunindextype neq, int nthreads, SUNContext sunctx) {
    return (nthreads > 1) ? N_VNew_OpenMP(neq, nthreads, sunctx)
                          : N_VNew_Serial(neq, sunctx);
}
} // namespace

// ---------------------------------------------------------------------------
// Construction / Destruction
// ---------------------------------------------------------------------------
FredIdaSolverBase::FredIdaSolverBase(const FuelRodGeometry& geom)
    : FredSolverBase(geom)
{}

FredIdaSolverBase::~FredIdaSolverBase() {
    if (m_LS)      SUNLinSolFree((SUNLinearSolver)m_LS);
    if (m_A)       SUNMatDestroy((SUNMatrix)m_A);
    if (m_ida_mem) IDAFree(&m_ida_mem);
}

// ---------------------------------------------------------------------------
// setupSundials — creates SUNDIALS context, vectors, IDA, linear solver,
//                 tolerances, id array, and algebraic suppression.
// ---------------------------------------------------------------------------
void FredIdaSolverBase::setupSundials(int max_steps, double t0) {
    const int neq = getSolverNeq();

    // Context
    SUNContext sunctx;
    if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != 0)
        throw std::runtime_error("FredIdaSolverBase: SUNContext_Create failed");
    m_sunctx = sunctx;

    // State vectors
    N_Vector y  = makeVector(neq, threads(), sunctx);
    N_Vector yp = makeVector(neq, threads(), sunctx);
    if (!y || !yp) throw std::runtime_error("FredIdaSolverBase: N_Vector allocation failed");
    m_y  = y;
    m_yp = yp;

    // Pack initial state into y; zero yp (IDACalcIC computes consistent yp).
    packSolverState(N_VGetArrayPointer(y));
    std::fill(N_VGetArrayPointer(yp), N_VGetArrayPointer(yp) + neq, 0.0);

    // IDA memory
    void* ida_mem = IDACreate(sunctx);
    if (!ida_mem) throw std::runtime_error("FredIdaSolverBase: IDACreate failed");
    m_ida_mem = ida_mem;

    IDAInit(ida_mem, getSolverResFn(), t0, y, yp);
    IDASetUserData(ida_mem, getSolverUserData());

    // Tolerances (scalar relative, per-variable absolute)
    N_Vector avtol = makeVector(neq, threads(), sunctx);
    N_VConst(m_atol, avtol);
    IDASVtolerances(ida_mem, m_rtol, avtol);
    N_VDestroy(avtol);

    // Dense direct linear solver
    SUNMatrix A = SUNDenseMatrix(neq, neq, sunctx);
    m_A = A;
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    m_LS = LS;
    IDASetLinearSolver(ida_mem, LS, A);

    if (max_steps > 0)
        IDASetMaxNumSteps(ida_mem, max_steps);
    if (m_hmax > 0)
        IDASetMaxStep(ida_mem, m_hmax);
    if (m_hinit > 0)
        IDASetInitStep(ida_mem, m_hinit);
    if (m_max_nonlin_iters > 0)
        IDASetMaxNonlinIters(ida_mem, m_max_nonlin_iters);
    if (m_max_ord > 0)
        IDASetMaxOrd(ida_mem, m_max_ord);

    // Mark differential (1) vs algebraic (0) variables.
    {
        N_Vector id_vec = makeVector(neq, threads(), sunctx);
        N_VConst(0.0, id_vec);
        fillSolverIdArray(N_VGetArrayPointer(id_vec), neq);
        IDASetId(ida_mem, id_vec);
        N_VDestroy(id_vec);
    }

    // Suppress algebraic variables from local truncation error test.
    IDASetSuppressAlg(ida_mem, 1);
}

// ---------------------------------------------------------------------------
// setupGapRoots
// ---------------------------------------------------------------------------
void FredIdaSolverBase::setupGapRoots(int nz, IDARootFn rootFn) {
    IDARootInit((void*)m_ida_mem, nz, rootFn);
}

// ---------------------------------------------------------------------------
// calcIC — IDACalcIC with get-consistent-IC on success, warning on failure.
// ---------------------------------------------------------------------------
bool FredIdaSolverBase::calcIC(double t_ic) {
    int rv = IDACalcIC((void*)m_ida_mem, IDA_YA_YDP_INIT, t_ic);
    if (rv < 0) {
        std::cerr << "FredIdaSolverBase: IDACalcIC failed (retval=" << rv
                  << ") — using user-supplied ICs\n";
        return false;
    }
    IDAGetConsistentIC((void*)m_ida_mem, (N_Vector)m_y, (N_Vector)m_yp);
    return true;
}

// ---------------------------------------------------------------------------
// reinitAfterEvent — re-initialise IDA after a gap open/close event.
//   Save y/yp, call IDAReInit, attempt IDACalcIC.
//   On CalcIC failure, restore saved state and reinit IDA with that.
// ---------------------------------------------------------------------------
void FredIdaSolverBase::reinitAfterEvent(double tret, double dtout) {
    const int neq = getSolverNeq();
    N_Vector y  = (N_Vector)m_y;
    N_Vector yp = (N_Vector)m_yp;

    std::vector<double> y_save(neq), yp_save(neq);
    std::copy(N_VGetArrayPointer(y),  N_VGetArrayPointer(y)  + neq, y_save.begin());
    std::copy(N_VGetArrayPointer(yp), N_VGetArrayPointer(yp) + neq, yp_save.begin());

    IDAReInit((void*)m_ida_mem, tret, y, yp);

    double dt_ic = (dtout > 0.0) ? std::min(dtout, 1.0) : 1.0;
    int rv = IDACalcIC((void*)m_ida_mem, IDA_YA_YDP_INIT, tret + dt_ic);
    if (rv >= 0) {
        IDAGetConsistentIC((void*)m_ida_mem, y, yp);
    } else {
        std::cerr << "FredIdaSolverBase: IDACalcIC after gap event (retval=" << rv
                  << ") — restoring pre-event state\n";
        std::copy(y_save.begin(),  y_save.end(),  N_VGetArrayPointer(y));
        std::copy(yp_save.begin(), yp_save.end(), N_VGetArrayPointer(yp));
        IDAReInit((void*)m_ida_mem, tret, y, yp);
    }
}

// ---------------------------------------------------------------------------
// Gap-event handlers — delegate to rodResiduals().
// ---------------------------------------------------------------------------
bool FredIdaSolverBase::isGapOpen(int j) {
    return rodResiduals().isGapOpen(j);
}

void FredIdaSolverBase::handleGapClosed(int j) {
    rodResiduals().setGapClosed(j,
        N_VGetArrayPointer((N_Vector)m_y),
        N_VGetArrayPointer((N_Vector)m_yp));
}

void FredIdaSolverBase::handleGapReopened(int j) {
    rodResiduals().setGapOpen(j, true);
}

// ---------------------------------------------------------------------------
// afterGapEventsHandled — default: re-initialise IDA after gap event.
// ---------------------------------------------------------------------------
void FredIdaSolverBase::afterGapEventsHandled(double tret, double dtout) {
    reinitAfterEvent(tret, dtout);
}

// ---------------------------------------------------------------------------
// runTimeLoop — shared IDA one-step loop with gap-event and output handling.
// ---------------------------------------------------------------------------
void FredIdaSolverBase::runTimeLoop(double tend, double dtout, bool all_steps, double t_start) {
    std::vector<int> rootsfound(m_geom.nz, 0);
    double t_next_out = t_start + dtout, tret = t_start;
    double t_loop_prev = t_start;  // tracks last accepted step time for afterAcceptedStep dt
    N_Vector y  = (N_Vector)m_y;
    N_Vector yp = (N_Vector)m_yp;

    // Snapshot housekeeping: copy timings so we can drain them as they are hit.
    std::vector<double> pending_snaps = m_snapshot_timings;
    m_snapshot_count = 0;

    while (true) {
        // Stop IDA exactly at the next output time so storeOutput is called
        // at the requested instant, not at whatever internal step IDA happens
        // to take (which can overshoot by a large factor in ONE_STEP mode).
        IDASetStopTime(m_ida_mem, t_next_out);
        int retval = IDASolve(m_ida_mem, t_next_out, &tret, y, yp, IDA_ONE_STEP);
        if (retval < 0) {
            std::cerr << "IDASolve failed at t=" << tret << " (retval=" << retval << ")\n";
            break;
        }

        // Physics hook: called at every accepted IDA step (subclasses override).
        // verbosity >= 5: time this hook so a wall-clock-heavy physics update
        // (e.g. GRSIS bubble evolution) can be distinguished from a legitimately
        // stiff IDA step-size collapse (verbosity >= 4 below) as the cause of a
        // slow run.
        const auto t_phys_start = std::chrono::steady_clock::now();
        afterAcceptedStep(tret, tret - t_loop_prev);
        if (m_verbosity >= 5) {
            const double phys_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t_phys_start).count();
            std::cout << "      [physics] afterAcceptedStep wall=" << phys_ms << " ms\n";
        }
        t_loop_prev = tret;

        // verbosity >= 4: per-accepted-step IDA solver diagnostics. Distinguishes
        // a legitimately stiff region (nsteps climbing fast, h shrinking, but
        // errTestFails/nonlinConvFails flat) from a solver actually struggling
        // (errTestFails/nonlinConvFails climbing) at every single internal step,
        // not just at dtout output points (logStepOutput only fires there).
        if (m_verbosity >= 4) {
            long nsteps = 0, nrevals = 0, netfails = 0, nncfails = 0;
            int    qlast = 0;
            double hlast = 0.0;
            IDAGetNumSteps(m_ida_mem, &nsteps);
            IDAGetNumResEvals(m_ida_mem, &nrevals);
            IDAGetNumErrTestFails(m_ida_mem, &netfails);
            IDAGetNumNonlinSolvConvFails(m_ida_mem, &nncfails);
            IDAGetLastOrder(m_ida_mem, &qlast);
            IDAGetLastStep(m_ida_mem, &hlast);
            std::cout << "    [ida] step#" << nsteps << " t=" << tret
                      << " h=" << hlast << " order=" << qlast
                      << " resEvals=" << nrevals
                      << " errTestFails=" << netfails
                      << " nonlinConvFails=" << nncfails << "\n";
        }

        bool is_root = (retval == IDA_ROOT_RETURN);
        if (is_root) {
            IDAGetRootInfo(m_ida_mem, rootsfound.data());
            bool event = false;
            for (int j = 0; j < m_geom.nz; ++j) {
                if (rootsfound[j] == -1 && isGapOpen(j)) {
                    std::cout << "  t=" << tret << " s: gap CLOSED at layer " << j << "\n";
                    handleGapClosed(j);
                    event = true;
                } else if (rootsfound[j] == 1 && !isGapOpen(j)) {
                    std::cout << "  t=" << tret << " s: gap REOPENED at layer " << j << "\n";
                    handleGapReopened(j);
                    event = true;
                }
            }
            if (event) afterGapEventsHandled(tret, dtout);
        }

        bool at_dtout = (retval == IDA_TSTOP_RETURN ||
                         tret >= t_next_out - 1.0e-12 * dtout);
        bool saved_snap_this_step = false;

        if (at_dtout) {
            unpackCurrentState();
            storeOutput(tret);
            if (m_verbosity >= 3) {
                double h_next = 0.0;
                IDAGetCurrentStep(m_ida_mem, &h_next);
                logStepOutput(tret, h_next);
            }

            // Checkpoint: overwrite single fault-recovery file.
            if (!m_ckpt_prefix.empty())
                saveCheckpoint(tret);

            // User-specified snapshot timings: save and drain matched entries.
            if (!m_snapshot_prefix.empty()) {
                for (auto it = pending_snaps.begin(); it != pending_snaps.end(); ) {
                    if (std::abs(tret - *it) < 0.5 * dtout) {
                        ++m_snapshot_count;
                        std::string fname = m_snapshot_prefix + "_frame"
                            + std::to_string(m_snapshot_count) + ".snapshot";
                        saveSnapshot(tret, fname);
                        saved_snap_this_step = true;
                        it = pending_snaps.erase(it);
                    } else {
                        ++it;
                    }
                }
            }

            t_next_out += dtout;
        } else if (all_steps && !is_root) {
            unpackCurrentState();
            storeOutput(tret);
        }

        if (tret >= tend - 1.0e-12 * dtout) {
            // Always save a final snapshot (skip if we just saved one at this time).
            if (!m_snapshot_prefix.empty() && !saved_snap_this_step) {
                ++m_snapshot_count;
                std::string fname = m_snapshot_prefix + "_frame"
                    + std::to_string(m_snapshot_count) + ".snapshot";
                saveSnapshot(tret, fname);
            }
            break;
        }
    }
}

// ---------------------------------------------------------------------------
// runHotStart — see FredIdaSolverBase.hpp doc comment.
// ---------------------------------------------------------------------------
void FredIdaSolverBase::runHotStart() {
    RodResiduals& res = rodResiduals();
    const bool prev_irr = res.isIrradiationOn();
    res.setEnableIrradiation(false);
    res.setBcTimeOverride(true, 0.0);

    if (m_verbosity >= 1)
        std::cout << "  hot start: marching to steady state at t=0 boundary "
                     "conditions (irradiation physics off)\n";

    const int neq = getSolverNeq();
    std::vector<double> id(neq);
    fillSolverIdArray(id.data(), neq);

    N_Vector y  = (N_Vector)m_y;
    N_Vector yp = (N_Vector)m_yp;
    std::vector<double> y_prev(N_VGetArrayPointer(y), N_VGetArrayPointer(y) + neq);

    const double tol       = 1.0e-4;  // max |dy| (K for temperatures) between pseudo-steps
    const int    max_iters = 60;
    double t_ps = 0.0;
    double dt   = (m_hinit > 0.0) ? m_hinit : 1.0;
    bool   converged = false;

    for (int iter = 0; iter < max_iters; ++iter) {
        double t_target = t_ps + dt;
        double tret = t_ps;
        int rv = IDASolve(m_ida_mem, t_target, &tret, y, yp, IDA_NORMAL);
        if (rv < 0) {
            std::cerr << "  hot start: IDASolve failed at pseudo-t=" << tret
                      << " (retval=" << rv << ") — stopping march, using best "
                         "available state\n";
            break;
        }
        t_ps = tret;

        double* yv = N_VGetArrayPointer(y);
        double maxdiff = 0.0;
        for (int k = 0; k < neq; ++k)
            if (id[k] > 0.5) maxdiff = std::max(maxdiff, std::fabs(yv[k] - y_prev[k]));
        std::copy(yv, yv + neq, y_prev.begin());

        if (m_verbosity >= 4)
            std::cout << "    [hot start] pseudo-t=" << t_ps
                      << " s  max|dy|=" << maxdiff << "\n";

        if (maxdiff < tol) { converged = true; break; }
        dt *= 2.0;
    }

    if (m_verbosity >= 1)
        std::cout << "  hot start: " << (converged ? "converged" : "reached iteration cap")
                  << " after " << t_ps << " s of pseudo-time\n";

    res.setBcTimeOverride(false);
    res.setEnableIrradiation(prev_irr);

    // Re-anchor IDA at t=0 with the converged hot state as the new initial
    // condition. IDASolve requires monotonically increasing tout, so the
    // pseudo-time march above cannot simply "become" the real timeline —
    // the caller's real time loop must start counting from t=0 again.
    IDAReInit((void*)m_ida_mem, 0.0, y, yp);
    calcIC(1.0);
}

// ---------------------------------------------------------------------------
// applyRestartToIDA — inject checkpoint state into a freshly initialised IDA.
// Call AFTER setupSundials(max_steps, t_start).
// ---------------------------------------------------------------------------
void FredIdaSolverBase::applyRestartToIDA(double t_start) {
    const int nz  = m_geom.nz;

    // Overwrite y/yp with checkpoint data
    double* y  = N_VGetArrayPointer((N_Vector)m_y);
    double* yp = N_VGetArrayPointer((N_Vector)m_yp);
    std::copy(m_restart_y .begin(), m_restart_y .end(), y);
    std::copy(m_restart_yp.begin(), m_restart_yp.end(), yp);

    // Restore gap state
    for (int j = 0; j < nz; ++j)
        rodResiduals().restoreGapState(j, m_restart_gap_open[j], m_restart_axial_offset[j]);

    // Re-initialise IDA with the injected state at the restart time.
    // Use t_start (caller's local copy) not m_restart_time which may already
    // have been reset to -1 in the run() method.
    IDAReInit((void*)m_ida_mem, t_start, (N_Vector)m_y, (N_Vector)m_yp);
}

} // namespace fred
