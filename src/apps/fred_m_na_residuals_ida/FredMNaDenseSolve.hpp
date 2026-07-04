#pragma once
// -----------------------------------------------------------------------
// FredMNaDenseSolve — small self-contained dense linear algebra + Newton
// helper for FredMNaSolver's one-step (fixed-dt, backward-Euler, "always
// accept") time integrator.
//
// FRED-M-Na no longer drives its time loop through SUNDIALS IDA (see
// FredMNaSolver::runOneStepLoop). Reusing IDA's own dense Newton/linear-
// algebra internals for the replacement would defeat the point, so this
// file provides a minimal from-scratch Gaussian-elimination-with-partial-
// pivoting solve and a generic dense Newton driver (finite-difference
// Jacobian), used to solve each axial layer's local nonlinear block
// (~100-150 unknowns) at every fixed time step. No SUNDIALS dependency.
// -----------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>

namespace fred {

// Solves A*x = b in place (A: row-major n*n, overwritten; b: overwritten
// with the solution). Returns false if a pivot is numerically singular
// (caller should just keep the current iterate — "always accept" semantics,
// no retry/shrink logic here).
//
// FRED-M-Na's residual rows mix wildly different physical units in a single
// dense block (e.g. a power-density AE row ~1e9 W/m3 next to a strain
// compatibility row ~1e0) — an ordinary partial-pivoting solve on the raw
// matrix loses rows to roundoff (accumulated cancellation drives an entire
// row to <1e-300 well before any genuine rank deficiency), which legacy
// never encounters because it solves each physical quantity by direct
// substitution/its own small same-units matrix, never one mixed-units
// system. Row equilibration (scale every row to unit max-norm before
// elimination) fixes this cheaply without changing the solution.
inline bool denseSolveInPlace(std::vector<double>& A, std::vector<double>& b, int n) {
    for (int row = 0; row < n; ++row) {
        double rmax = 0.0;
        for (int k = 0; k < n; ++k) rmax = std::max(rmax, std::fabs(A[row*n + k]));
        if (rmax < 1.0e-300) return false; // row identically zero: truly singular
        const double s = 1.0 / rmax;
        for (int k = 0; k < n; ++k) A[row*n + k] *= s;
        b[row] *= s;
    }

    for (int col = 0; col < n; ++col) {
        int piv = col;
        double best = std::fabs(A[col*n + col]);
        for (int row = col+1; row < n; ++row) {
            double v = std::fabs(A[row*n + col]);
            if (v > best) { best = v; piv = row; }
        }
        if (best < 1.0e-12) return false;
        if (piv != col) {
            for (int k = 0; k < n; ++k) std::swap(A[col*n + k], A[piv*n + k]);
            std::swap(b[col], b[piv]);
        }
        const double diag = A[col*n + col];
        for (int row = col+1; row < n; ++row) {
            const double factor = A[row*n + col] / diag;
            if (factor == 0.0) continue;
            for (int k = col; k < n; ++k) A[row*n + k] -= factor * A[col*n + k];
            b[row] -= factor * b[col];
        }
    }
    for (int row = n-1; row >= 0; --row) {
        double sum = b[row];
        for (int k = row+1; k < n; ++k) sum -= A[row*n + k] * b[k];
        const double diag = A[row*n + row];
        if (std::fabs(diag) < 1.0e-300) return false;
        b[row] = sum / diag;
    }
    return true;
}

// Generic dense Newton solve: finds y such that residualFn(y) ~= 0, via
// finite-difference Jacobian + denseSolveInPlace. Always returns whatever
// iterate it ends on (converged or not) after max_iter — no step-size
// control, no rejection — matching the legacy FRED-M "always accept the
// current iterate" philosophy this integrator replicates. y is updated
// in place. Returns the infinity-norm of the last Newton correction (0 if
// it converged and took no further step), useful as an outer-loop
// convergence signal for the caller.
//
// scale: optional (nullable) per-unknown "typical physical magnitude" array
// (size n) — used for two things, both needed because this system mixes
// wildly different physical units in one dense block (e.g. a power-density
// unknown ~1e9 W/m3 next to a strain unknown ~1e-3):
//   1. Floors the finite-difference step. A plain relative step
//      (eps*|value|) silently underflows to a bitwise-zero derivative
//      whenever a large-scale unknown sits at exactly 0 (e.g. qqv before
//      power ramps up — perturbing 0 by a tiny epsilon is completely lost
//      against the ~1e9-scale residual it feeds, since that's under one
//      ULP at that magnitude).
//   2. Floors the per-component Newton-step damping limit (see below).
// Pass nullptr to fall back to a plain relative step/limit (fine for
// reasonably-uniformly-scaled systems; FRED-M-Na's mixed-units blocks need
// the real thing).
inline double denseNewtonSolve(int n,
                                const std::function<void(const double* y, double* r)>& residualFn,
                                double* y,
                                int max_iter,
                                double atol,
                                const double* scale = nullptr)
{
    std::vector<double> r(n), r_pert(n), J(static_cast<size_t>(n) * n), rhs(n);
    std::vector<double> y_trial(y, y + n);

    residualFn(y_trial.data(), r.data());
    double last_delta = 0.0;

    for (int it = 0; it < max_iter; ++it) {
        double rnorm = 0.0;
        for (int i = 0; i < n; ++i) rnorm = std::max(rnorm, std::fabs(r[i]));
        if (rnorm < atol) break;

        for (int k = 0; k < n; ++k) {
            const double save = y_trial[k];
            const double typical = (scale && scale[k] > 0.0) ? scale[k] : 1.0;
            const double h = std::max(1.0e-6 * typical, 1.0e-6 * std::fabs(save));
            y_trial[k] = save + h;
            residualFn(y_trial.data(), r_pert.data());
            y_trial[k] = save;
            for (int i = 0; i < n; ++i)
                J[static_cast<size_t>(i)*n + k] = (r_pert[i] - r[i]) / h;
        }
        for (int i = 0; i < n; ++i) rhs[i] = -r[i];

        if (!denseSolveInPlace(J, rhs, n))
            break; // singular: stop refining, accept current

        // Damp the raw Newton correction: this system is strongly nonlinear
        // (e.g. temperature-dependent conductivity feeding straight back
        // into the thermal residual it's being solved for), and the first
        // iteration from a cold, uniform-temperature start under full power
        // is a large excursion — an undamped full Newton step can overshoot
        // into an unphysical region (e.g. T runs to 1e7+ K) that the next
        // iteration's Jacobian can't recover from.
        //
        // Scale the WHOLE correction by a single shared factor (not an
        // independent per-component clip): clipping components
        // independently breaks the Newton direction (it stops being a
        // solution of the linearised system at all) and was found to
        // diverge outright on this problem's tightly-coupled temperature
        // chain. A single scalar factor keeps the direction intact and
        // only shortens the step length — standard damped-Newton/line-
        // search practice. Generous per-component limits (an order of
        // magnitude above the variable's own current-or-typical size) still
        // stop genuine blow-ups without starving legitimately large first-
        // iteration moves (e.g. qqv jumping from 0 to ~1e9 under a power
        // step). This does not change what "always accept" means (the
        // OUTER fixed-step is still always accepted whatever it converges
        // to) — it just keeps the INNER Newton iteration well-behaved.
        double damp = 1.0;
        for (int k = 0; k < n; ++k) {
            const double limit = 10.0 * std::max(std::fabs(y_trial[k]),
                                                  scale ? scale[k] : 1.0);
            if (limit > 0.0 && std::fabs(rhs[k]) > limit)
                damp = std::min(damp, limit / std::fabs(rhs[k]));
        }

        last_delta = 0.0;
        for (int k = 0; k < n; ++k) {
            const double step = damp * rhs[k];
            y_trial[k] += step;
            last_delta = std::max(last_delta, std::fabs(step));
        }
        residualFn(y_trial.data(), r.data());
    }

    std::copy(y_trial.begin(), y_trial.end(), y);
    return last_delta;
}

} // namespace fred
