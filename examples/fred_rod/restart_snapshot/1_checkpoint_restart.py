"""
FRED-ROD restart_snapshot — Script 1: Restart from checkpoint.

Scenario: the simulation in script 0 ran to t_half and then crashed.
The checkpoint file (rod_checkpoint.chk) contains the state at t_half.
This script loads that checkpoint and continues the simulation to tend.

Simulation time is preserved (restart from t_half, not reset to 0).

Run script 0 first to generate the checkpoint file.
"""

import sys, os
import numpy as np

BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import fred_rod as fred

WORK_DIR = os.path.join(os.path.dirname(__file__), 'ckpt')

# ── Simulation parameters (must match script 0) ───────────────────────────────
tend   = 5000.0
dtout  = 1000.0
t_half = 3000.0

# For this demonstration we re-run from 0 to t_half to generate the checkpoint.
def make_solver():
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
    s.set_power_density_history([0.0, tend], [2.0e8, 2.0e8])
    s.set_coolant_temperature(  [0.0, tend], [600.0, 600.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(600.0)
    s.set_initial_gas_pressure(0.1)
    s.set_tolerances(1.0e-5, 1.0e-7)
    return s


ckpt_prefix = os.path.join(WORK_DIR, "rod")
ckpt_file   = os.path.join(WORK_DIR, "rod_checkpoint.chk")

# ── Step A: run 0 → t_half (simulating the partial run before crash) ──────────
print("=" * 60)
print(f"Step A: partial run 0 → {t_half:.0f} s  (checkpoints enabled)")
print("=" * 60)

s_partial = make_solver()
s_partial.set_checkpoint_prefix(ckpt_prefix)
s_partial.run(t_half, dtout)

print(f"Partial run done.  Checkpoint written: {os.path.basename(ckpt_file)}")

# ── Step B: restart from checkpoint at t_half, continue to tend ───────────────
print("\n" + "=" * 60)
print(f"Step B: load checkpoint → restart from t_half to tend")
print("=" * 60)

s_restart = make_solver()
s_restart.set_checkpoint_prefix(ckpt_prefix)   # continue checkpointing
s_restart.load_checkpoint(ckpt_file)
print(f"Loaded checkpoint.  Restart time: {s_restart.restart_time():.0f} s")

s_restart.run(tend, dtout)

t_restart   = np.array(s_restart.time_points())
T_restart   = np.array(s_restart.peak_fuel_temperature())
gap_restart = np.array(s_restart.gap_width())
rco_restart = np.array(s_restart.clad_outer_radius())

print(f"\nRestart run: {len(t_restart)} output points, "
      f"t=[{t_restart[0]:.0f}, {t_restart[-1]:.0f}] s")

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
    print("Verification: checkpoint restart vs reference (script 0)")
    print("=" * 60)

    tol_t = 1.0
    matches = [(i, j, t) for i, t in enumerate(t_ref)
               for j, tr in enumerate(t_restart) if abs(t - tr) < tol_t]
    print(f"Overlapping time points: {len(matches)}")

    def check(name, a_ref, a_rst, rtol=1e-6):
        errs = [abs(a_ref[i] - a_rst[j]) / max(abs(a_ref[i]), abs(a_rst[j]), 1e-30)
                for i, j, _ in matches]
        mx = max(errs) if errs else 0.0
        ok = mx < rtol
        print(f"  {name:<35} max_rel_err={mx:.2e}  [{'PASS' if ok else 'FAIL'}]")
        return ok

    ok = True
    ok &= check("peak_fuel_temperature [K]", T_ref,   T_restart)
    ok &= check("gap_width [m]",             gap_ref, gap_restart)
    ok &= check("clad_outer_radius [m]",     rco_ref, rco_restart)

    print()
    print("OVERALL:", "PASS" if ok else "FAIL",
          "— checkpoint restart reproduces the reference exactly." if ok
          else "— mismatch detected.")
