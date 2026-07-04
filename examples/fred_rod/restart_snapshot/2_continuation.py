"""
FRED-ROD restart_snapshot — Script 2: Continuation via snapshot.

Scenario: "irradiation" phase (0 → t_half) followed by a "transient" phase.

  Phase 1: run from t0 to t_half.  A snapshot is automatically saved at the
           end of the run (rod_frame1.snapshot).  This represents the fuel
           state at the end of the irradiation.

  Phase 2: load the snapshot.  Simulation time resets to 0.0 (fresh start)
           while the physical fuel state (temperature, strain, gap) is preserved.
           Run from 0 to (tend - t_half).

Combined output:
  Phase 1 times:  [0, dtout, 2*dtout, ..., t_half]
  Phase 2 times:  [0, dtout, ..., tend-t_half]  → offset by t_half for plotting

The combined trajectory should match the reference run (script 0).
"""

import sys, os
import numpy as np

BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import fred_rod as fred

WORK_DIR = os.path.join(os.path.dirname(__file__), 'ckpt')
os.makedirs(WORK_DIR, exist_ok=True)

# ── Simulation parameters ─────────────────────────────────────────────────────
tend   = 5000.0
dtout  = 1000.0
t_half = 3000.0   # end of phase 1 / start of phase 2
t_phase2 = tend - t_half   # duration of phase 2


def make_solver(t_start=0.0, t_end=tend):
    """Build solver with boundary conditions covering [t_start, t_end]."""
    g = fred.FuelRodGeometry()
    g.nf  = 4;  g.nc  = 3;  g.nz  = 2
    g.rfi0 = 0.0;  g.rfo0 = 4.0e-3
    g.rci0 = 4.1e-3;  g.rco0 = 4.7e-3
    g.dz0  = [0.15, 0.15]
    g.vgp  = 1.0e-6;  g.ruff = 5e-6;  g.rufc = 5e-6
    g.build()

    fuel = fred.UO2(reference_density=10400.0)
    clad = fred.AIM1(reference_density=7900.0)
    gap  = fred.He()

    s = fred.FredRodSolver(g, fuel, clad, gap)
    # Boundary conditions are the same in both phases.
    t_bc = [0.0, max(t_end, tend)]
    s.set_power_density_history(t_bc, [2.0e8, 2.0e8])
    s.set_coolant_temperature(  t_bc, [600.0, 600.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(600.0)
    s.set_initial_gas_pressure(0.1)
    s.set_tolerances(1.0e-5, 1.0e-7)
    return s


snap_prefix = os.path.join(WORK_DIR, "rod")
snap_file   = os.path.join(WORK_DIR, "rod_frame1.snapshot")  # auto-saved at t_half

# ── Phase 1: run 0 → t_half, snapshot saved at end ───────────────────────────
print("=" * 60)
print(f"Phase 1: run 0 → {t_half:.0f} s  (snapshot saved at end)")
print("=" * 60)

s1 = make_solver()
# No user-specified snapshot times; the final state is always saved.
s1.set_snapshot_prefix(snap_prefix)
s1.run(t_half, dtout)

t1   = np.array(s1.time_points())
T1   = np.array(s1.peak_fuel_temperature())
gap1 = np.array(s1.gap_width())
rco1 = np.array(s1.clad_outer_radius())

print(f"Phase 1 done: {len(t1)} points, t=[{t1[0]:.0f}, {t1[-1]:.0f}] s")
print(f"  Snapshot saved: {os.path.basename(snap_file)}")

# ── Phase 2: load snapshot, run 0 → t_phase2 ─────────────────────────────────
print("\n" + "=" * 60)
print(f"Phase 2: load snapshot → run 0 → {t_phase2:.0f} s (time reset)")
print("=" * 60)

s2 = make_solver(t_end=t_phase2)
s2.load_snapshot(snap_file)
print(f"  restart_time = {s2.restart_time():.1f} s  (0 = time reset confirmed)")

s2.run(t_phase2, dtout)

t2   = np.array(s2.time_points())
T2   = np.array(s2.peak_fuel_temperature())
gap2 = np.array(s2.gap_width())
rco2 = np.array(s2.clad_outer_radius())

print(f"Phase 2 done: {len(t2)} points, t=[{t2[0]:.0f}, {t2[-1]:.0f}] s")

# ── Combined trajectory (offset phase-2 times by t_half) ─────────────────────
t_combined   = np.concatenate([t1, t2[1:] + t_half])
T_combined   = np.concatenate([T1, T2[1:]])
gap_combined = np.concatenate([gap1, gap2[1:]])
rco_combined = np.concatenate([rco1, rco2[1:]])

print(f"\nCombined: {len(t_combined)} points, "
      f"t=[{t_combined[0]:.0f}, {t_combined[-1]:.0f}] s")

# ── Verify against reference ──────────────────────────────────────────────────
ref_file = os.path.join(WORK_DIR, "ref_times.npy")
if not os.path.exists(ref_file):
    print("\n[INFO] Reference data not found — run script 0 first for comparison.")
else:
    t_ref   = np.load(os.path.join(WORK_DIR, "ref_times.npy"))
    T_ref   = np.load(os.path.join(WORK_DIR, "ref_T.npy"))
    gap_ref = np.load(os.path.join(WORK_DIR, "ref_gap.npy"))
    rco_ref = np.load(os.path.join(WORK_DIR, "ref_rco.npy"))

    print("\n" + "=" * 60)
    print("Verification: continuation vs reference (script 0)")
    print("=" * 60)

    tol_t = 1.0
    matches = [(i, j, t) for i, t in enumerate(t_ref)
               for j, tc in enumerate(t_combined) if abs(t - tc) < tol_t]
    print(f"Overlapping time points: {len(matches)}")

    def check(name, a_ref, a_comb, rtol=1e-6):
        errs = [abs(a_ref[i] - a_comb[j]) / max(abs(a_ref[i]), abs(a_comb[j]), 1e-30)
                for i, j, _ in matches]
        mx = max(errs) if errs else 0.0
        ok = mx < rtol
        print(f"  {name:<35} max_rel_err={mx:.2e}  [{'PASS' if ok else 'FAIL'}]")
        return ok

    ok = True
    ok &= check("peak_fuel_temperature [K]", T_ref,   T_combined)
    ok &= check("gap_width [m]",             gap_ref, gap_combined)
    ok &= check("clad_outer_radius [m]",     rco_ref, rco_combined)

    print()
    print("OVERALL:", "PASS" if ok else "FAIL",
          "— continuation matches reference exactly." if ok
          else "— mismatch detected.")
