# FRED 2.0 Platform Architecture

## Source tree layout

```
src/
  platform/                  Abstract base classes and shared solvers
    FuelPelletMaterial.hpp   Abstract fuel pellet interface
    CladdingMaterial.hpp     Abstract cladding interface
    GapMaterial.hpp          Abstract gap conductance interface
    GapModel.hpp / .cpp      Free functions: gas conductance, contact conductance
    HeatConduction.hpp / .cpp  FVM thermal solver (cylindrical, per axial layer)
    StressStrain.hpp / .cpp    Thermo-elastic mechanics solver
    FuelRodGeometry.hpp      Geometry struct + build()
    FuelRodState.hpp         State struct (temperatures, strains, stresses)
    Constants.hpp            Physical constants

  apps/
    fred_rod/                FRED-ROD: pure thermo-mechanical, no irradiation
      fuelpelletmaterial/
        UO2.hpp / .cpp       UO2 at zero burnup, as-fabricated density
        DummyFuelPellet.hpp / .cpp  Constant-property benchmark material
      claddingmaterial/
        AIM1.hpp / .cpp      AIM1 austenitic stainless steel
        DummyCladding.hpp / .cpp    Constant-property benchmark material
      gapmaterial/
        DummyGapMaterial.hpp Linear h(T) = 500 + 2.5*T [W/(m²·K)]
      FredRodResiduals.hpp / .cpp  DAE residual assembly
      FredRodSolver.hpp / .cpp     High-level solver (wraps SUNDIALS IDA)
      bindings.cpp                 Python/pybind11 module

    fred_ox/                 FRED-OX: mixed oxide with irradiation (in development)
      FredOxGapMaterial.hpp / .cpp  Gap conductance with fission gas mixture
      FredOxResiduals.hpp           Skeletal residual class (TODO: implement)
```

## Material inheritance hierarchy

### Fuel pellet

```
FuelPelletMaterial  (src/platform/)
  — thermalConductivity(T)         temperature only; no burnup/density
  — heatCapacity(T)
  — thermalExpansionStrain(T)
  — youngsModulus(T, density)      density kept: porosity correction
  — poissonRatio()
  — referenceDensity()
  — theoreticalDensity()
  — meltingTemperature()

  UO2  (apps/fred_rod/fuelpelletmaterial/)
    — Philipponneau (1992) at zero burnup, stored as-fabricated density
    — Valid for fresh fuel / out-of-pile tests (FRED-ROD use case)

  DummyFuelPellet  (apps/fred_rod/fuelpelletmaterial/)
    — All properties constant; for benchmarks and unit tests

  [FRED-OX] FredOxFuelPelletMaterial  (apps/fred_ox/) — TO BE IMPLEMENTED
    — Inherits from UO2 (or directly from FuelPelletMaterial)
    — Overrides thermalConductivity(T) to read current burnup/density from
      solver state, enabling the full irradiation-dependent correlation
    — Adds swelling strain, creep strain, fission gas production
```

### Cladding

```
CladdingMaterial  (src/platform/)
  — thermalConductivity(T)
  — heatCapacity(T)
  — thermalExpansionStrain(T)
  — youngsModulus(T)
  — poissonRatio()
  — meyerHardness(T)
  — referenceDensity()

  AIM1  (apps/fred_rod/claddingmaterial/)
    — AIM1 austenitic stainless steel; KIT/CEA correlations

  DummyCladding  (apps/fred_rod/claddingmaterial/)
    — Constant properties; for benchmarks

  [FRED-OX] FredOxCladdingMaterial  (apps/fred_ox/) — TO BE IMPLEMENTED
    — Inherits from AIM1 (or CladdingMaterial directly)
    — Adds irradiation creep, void swelling
```

### Gap conductance

```
GapMaterial  (src/platform/)
  — gapConductance(T)   open-gap h [W/(m²·K)] as a function of temperature only

  DummyGapMaterial  (apps/fred_rod/gapmaterial/)
    — h(T) = 500 + 2.5*T  [W/(m²·K)]
    — Calibrated to He-filled gap of ~100 µm (Lanning-Hann model)
    — Default for FRED-ROD

  FredOxGapMaterial  (apps/fred_ox/)
    — Inherits from GapMaterial
    — Holds fission gas state: m_mu0 (He fill gas), m_fgrel (released FG),
      m_gpres, m_burnup
    — gapConductance(T) calls gasMixtureConductivity() for k_gas from the
      Xe/Kr/He fission gas mixture, then applies the Lanning-Hann model
    — Solver updates state via setFissionGasRelease() / setBurnup() at each
      output step before the residual is evaluated
```

## Design rules for adding FRED-OX materials

FRED-OX materials must inherit from the corresponding FRED-ROD concrete class
(or from the platform abstract class) and override only the methods that change
under irradiation.  The guiding principle is that FRED-OX is a superset of
FRED-ROD: every property that is valid for fresh fuel at as-fabricated density
is inherited unchanged.

Checklist when implementing a new FRED-OX material:

1. Inherit from the FRED-ROD class (e.g. `UO2`, `AIM1`) so correlations for
   thermal, elastic, and expansion properties are available by default.
2. Add irradiation state as private members (burnup, fast fluence, fgrel, …).
3. Override only the methods that depend on that state.
4. Provide setter methods (`setBurnup`, `setFluence`, …) that the solver calls
   at each time step before assembling the residual.
5. Do not change the abstract interface in `src/platform/` — new irradiation
   parameters belong in the FRED-OX subclass, not in the base.

## Gap conductance split: open vs. closed gap

`HeatConduction` uses two code paths for the gap:

- **Open gap** (`s.gapOpen == true`):  
  `h = m_gap_mat.gapConductance(T_avg)`  
  The entire gas + radiation model is encapsulated in `GapMaterial`. FRED-ROD
  uses `DummyGapMaterial`; FRED-OX uses `FredOxGapMaterial`.

- **Closed gap** (`s.gapOpen == false`):  
  `h = computeContactConductance(...)` from `GapModel.hpp`.  
  This is a platform-level free function (Ross & Stoute model) that is shared
  across all apps and does not depend on fission gas state.

### FRED-M-Na: three-state (open/soft/clos) gap-contact model

FRED-ROD and FRED-OX use a **two-state** gap model (`AxialLayerState::gapOpen`,
a bool), with transitions detected by a SUNDIALS IDA root (`RodResiduals::gapRoot`,
`FredSolverBase::setupGapRoots`). This is unchanged and remains the shared
platform mechanism — see `src/platform/GapStateManager.{hpp,cpp}` and
`src/platform/RodResiduals.{hpp,cpp}`.

FRED-M-Na (U-Pu-Zr metallic fuel in sodium) instead implements legacy
FRED-M's **three-state** `open`/`soft`/`clos` gap-contact model, because
anisotropic axial swelling produces localized fuel-clad contact (`soft`) well
before the bulk annular gap geometrically closes (`clos`) — a real physical
regime with no two-state analogue. This model is entirely app-specific and
does **not** touch any platform file:

- **State tracking**: `FredMNaLayerState::flag` (`"open"`/`"soft"`/`"clos"`,
  `apps/fred_m_na/FredMNaState.hpp`), a monotonic ratchet (never reopens for
  UPuZr fuel, matching legacy `Baseir.for`).
- **Detection**: no IDA root-finder is registered for FRED-M-Na
  (`setupGapRoots` is not called). Instead, `FredMNaSolver::afterAcceptedStep`
  checks gap geometry directly on every accepted backward-Euler step (the
  same cadence legacy's fixed-step driver uses) and drives the ratchet via
  `FredMNaGapBehavior::update` (`apps/fred_m_na/FredMNaGapBehavior.{hpp,cpp}`).
- **Mechanics**: `apps/fred_m_na/FredMNaStressStrain.{hpp,cpp}` is a
  `.cpp`-only fork of `platform/StressStrain.cpp` (composition, not
  subclassing — `StressStrain::computeResiduals` is non-virtual). It reads
  `FredMNaLayerState::flag` instead of `AxialLayerState::gapOpen` at the 4 BC
  branch sites (fuel axial force balance, fuel outer radial BC, cladding
  axial/interface BC, cladding inner radial BC — `soft` groups with `clos` at
  the first three, but with `open` at the cladding inner radial BC, since
  contact is localized rather than full annular), and reads directional
  swelling accumulators `efsz`/`efsh`/`efsr` (also `FredMNaLayerState`-only)
  instead of the isotropic `efs/3` term.
- **Residual wiring**: `FredMNaResiduals::computeResidualsMNa` /
  `idaResidualMNa` reimplement `RodResiduals::computeResiduals`'s per-layer
  loop (which is non-virtual and platform-shared) to call
  `FredMNaStressStrain` instead of the inherited `StressStrain m_mech`.

`src/platform/FuelRodState.hpp`, `StressStrain.{hpp,cpp}`, `GapStateManager.{hpp,cpp}`,
and `RodResiduals.{hpp,cpp}` are unmodified by any of the above — FRED-ROD and
FRED-OX numerical output is unaffected by construction (verified bit-identical
against pre-refactor HDF5 baselines).
### FRED-M-Na: one-step backward-Euler integrator (no SUNDIALS IDA)

FRED-ROD and FRED-OX advance in time via `FredSolverBase::runTimeLoop`,
which drives SUNDIALS IDA (`IDASolve(..., IDA_ONE_STEP)`) — an adaptive-
order BDF corrector that shrinks the step and retries on a non-converged
Newton iteration. Legacy FRED-M (`Baseir.for`) does not do this: it takes a
**fixed** step `dt = dtout` and iterates a Picard (successive-substitution)
outer sweep, capped at a fixed iteration count (`niter > 1000`), and
**always accepts whatever iterate it ends on** — there is no step-size
shrink-and-retry path in the legacy code at all (`Baseir.for:367-371`,
`goto 100` on the cap being hit simply continues with the current state).
This mismatch was the root cause of repeated `IDASolve` corrector-
convergence-failure crashes in FRED-M-Na (reinit-event clustering around
the `open`/`soft`/`clos` gap-contact transition, ~t=94-165 d in the 200-day
benchmark) that could not be fixed by loosening IDA's `rtol`/`atol`, since
those only change which *converged* steps are accepted, not whether a
non-converged Newton iterate gets force-accepted.

FRED-M-Na therefore does **not** call `IDACreate`/`IDAInit`/`IDASolve`/
`IDAReInit` at all. `FredMNaSolver::run()` still allocates a SUNDIALS
`N_Vector` for `m_y`/`m_yp` (inherited from `FredSolverBase`) purely as a
plain data buffer, so the existing checkpoint/snapshot/HDF5-restart
plumbing keeps working unmodified — but the actual time-stepping is its
own fixed-step, backward-Euler, "always accept" driver, entirely inside
`src/apps/fred_m_na/`:

- **`FredMNaSolver::runOneStepLoop`** (replaces `FredSolverBase::runTimeLoop`
  for this app): fixed `dt` (defaults to `dtout`, or `set_step_size(...)`),
  calling `solveStepBackwardEuler` each step, then the same
  `afterAcceptedStep`/`storeOutput`/`logStepOutput`/checkpoint/snapshot
  hooks `runTimeLoop` used to call — so all existing per-step physics
  (GRSIS, Zr redistribution, cladding wastage, the open/soft/clos ratchet)
  is unchanged, just driven by this loop instead of IDA. There is no
  root-finding/`reinitAfterEvent` dispatch: with no adaptive step/order
  history to invalidate, gap-state transitions applied in
  `afterAcceptedStep` are simply picked up by the next fixed step's solve.
- **`FredMNaSolver::solveStepBackwardEuler`**: per step, an outer Picard
  sweep (capped at `MAX_OUTER=60`, in practice converging in 2-4
  iterations) — the 3 global rows (`fggen`/`fgrel`/`gpres`) are updated via
  `FredMNaResiduals::computeGlobalUpdate` (direct substitution; `gpres` is
  additionally solved jointly with each layer for tighter coupling, since
  its defining equation depends on each layer's own deformed geometry), and
  each axial layer's local block is solved independently via
  `newtonSolveLayer` (layers only couple through the 3 rod-scalar globals —
  see `FredMNaResiduals::computeLayerResiduals`). The sweep is **always
  accepted** whatever it ends on after `MAX_OUTER` iterations — no dt
  shrink, no retry — matching legacy's `goto 100` escape-and-continue.
- **`FredMNaDenseSolve.hpp`**: a small, self-contained, SUNDIALS-free
  Gaussian-elimination-with-row-equilibration linear solve
  (`denseSolveInPlace`) and a generic damped finite-difference Newton
  driver (`denseNewtonSolve`) used by `newtonSolveLayer` to solve each
  layer's local nonlinear block (~100-150 unknowns). Row equilibration and
  per-unknown finite-difference-step/damping scales
  (`FredMNaResiduals::fillFdScale`) are needed because a single layer's
  block mixes wildly different physical units (e.g. ~1e9 W/m³
  power-density rows next to ~1e-3 strain rows) in one dense system.

`src/platform/FredSolverBase.{hpp,cpp}` is unmodified by this — FRED-ROD
and FRED-OX still drive through `runTimeLoop`/IDA exactly as before
(verified bit-identical against pre-refactor HDF5 baselines). `FredMNaResiduals::computeResidualsMNa`/`idaResidualMNa` (the old full-vector
IDA residual entry point used before this change) still compile and remain
available for offline debugging/Jacobian sanity checks, but are no longer
on the hot path.

Net effect on the reference benchmark
(`examples/fred_m_na/timpano_SAS_benchmark/base_irradiation`, 2176 days):
what previously either crashed (`IDASolve` corrector failure around
t=163-166 d in the 200-day case) or took 4+ hours (full 2176-day case) now
completes in ~1.5 s (200 d) / ~80 s (2176 d), with end-of-run agreement
against legacy FRED-M of `xwast` 45.2 vs 45.8 µm and FGR 94.7% vs 93.6%.
Plenum gas pressure remains a known, unresolved discrepancy (4.6 vs
6.6 MPa at t=2176 d) — pre-existing, not introduced by this change; see
`scripts/session_state.md` for the current investigation status.
