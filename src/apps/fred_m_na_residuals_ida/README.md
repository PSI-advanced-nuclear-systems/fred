# fred_m_na_residuals_ida — reference snapshot (do not build)

This directory is a **reference-only copy** of `src/apps/fred_m_na` taken on
2026-07-02, immediately before the "remove dead SUNDIALS/IDA scaffolding"
refactor (see `scripts/session_state.md`, task "remove all sundials based
variables ... for fred-m-na which does not use this").

## Why this snapshot exists

By this point FRED-M-Na had already been rewritten to use its own fixed-step
backward-Euler "one-step" integrator (`FredMNaSolver::runOneStepLoop` /
`solveStepBackwardEuler` / `newtonSolveLayer`) instead of SUNDIALS IDA. But
the app still carried a full, unused SUNDIALS/IDA integration surface,
required only because `FredSolverBase` (the shared platform base class)
declared several pure-virtual "solver setup hooks" —
`getSolverResFn()` / `getSolverUserData()` / `fillSolverIdArray()` /
`packSolverState()` — that assumed **every** application drives itself
through `FredSolverBase::runTimeLoop()` / `setupSundials()`. FRED-M-Na never
calls either of those, so:

- `FredMNaResiduals::computeResidualsMNa` / `idaResidualMNa` (a full
  SUNDIALS `IDAResFn`-shaped residual entry point) had zero real callers —
  it existed solely so `FredMNaSolver::getSolverResFn()` had something to
  return to satisfy the base class's pure-virtual contract.
- `FredMNaSolver::handleGapClosed` / `handleGapReopened` overrides were
  likewise dead: both are only ever invoked from
  `FredSolverBase::runTimeLoop()`'s IDA root-event dispatch, which FRED-M-Na
  never reaches (no root-finder is registered for this app —
  `handleGapReopened`'s own body was `assert(false)`, documented as
  unreachable).
- `set_tolerances` / `set_max_step` / `set_init_step` /
  `set_max_nonlin_iters` / `set_max_order` were exposed on the Python
  binding for API-compatibility only; their docstrings already said
  "Unused by FRED-M-Na (no SUNDIALS IDA integrator) ... has no effect".

This snapshot preserves that code **as it was**, for reference, before the
follow-up refactor:
1. split `FredSolverBase` into a slim base (checkpoint/HDF5/output
   plumbing genuinely shared by every app) plus a new `FredIdaSolverBase`
   (the SUNDIALS-only pieces), so only FRED-ROD/FRED-OX — the two apps that
   actually call SUNDIALS IDA — are required to implement the IDA-shaped
   virtual hooks;
2. deleted the now-provably-dead `computeResidualsMNa` / `idaResidualMNa`,
   the dead `handleGapClosed`/`handleGapReopened` overrides, and the
   no-op tolerance/step-size Python bindings from the live
   `src/apps/fred_m_na/`.

## Status

**Not part of the build.** This directory is intentionally omitted from the
Makefile's source lists and is not compiled by `make all`. It is kept purely
so a future session can diff the pre-refactor design against the live code
if a question comes up about why something changed.
