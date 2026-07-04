"""
FRED-OX restart_snapshot — Script 1: Restart from checkpoint.

Scenario: simulation crashed after day 1.  The checkpoint (ox_checkpoint.chk)
contains the fuel state at the end of the first-day run.  This script loads
the checkpoint and resumes to tend.  Simulation time is preserved (not reset).

Run script 0 first to generate the reference data (optional, for verification).
"""

import sys, os
import numpy as np

BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import fred_ox as fred

WORK_DIR = os.path.join(os.path.dirname(__file__), 'ckpt')
os.makedirs(WORK_DIR, exist_ok=True)

tend   = 3.0 * 86400
dtout  = 86400.0
t_half = 86400.0

ckpt_prefix = os.path.join(WORK_DIR, "ox")
ckpt_file   = os.path.join(WORK_DIR, "ox_checkpoint.chk")


def make_solver():
    g = fred.FuelRodGeometry()
    g.nf  = 4;  g.nc  = 3;  g.nz  = 1
    g.rfi0 = 0.0;  g.rfo0 = 4.0e-3
    g.rci0 = 4.1e-3;  g.rco0 = 4.7e-3
    g.dz0  = [0.20]
    g.vgp  = 1.5e-6;  g.ruff = 5e-6;  g.rufc = 5e-6
    g.build()

    mox  = fred.FredOxMOX(pu_content=0.20, rof0=10200.0, sto0=1.97)
    clad = fred.FredOxT91(reference_density=7750.0)
    gap  = fred.FredOxGapMaterial()

    s = fred.FredOxSolver(g, mox, clad, gap)
    s.set_power_density_history([0.0, tend], [1.5e8, 1.5e8])
    s.set_coolant_temperature(  [0.0, tend], [620.0, 620.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(620.0)
    s.set_initial_gas_pressure(0.1)
    s.set_tolerances(1.0e-5, 1.0e-7)
    return s


# ── Step A: run 0 → t_half (partial run, simulates pre-crash) ────────────────
print("=" * 60)
print(f"Step A: partial run 0 → day 1 with checkpoints")
print("=" * 60)

s_partial = make_solver()
s_partial.set_checkpoint_prefix(ckpt_prefix)
s_partial.run(t_half, dtout)
print(f"Partial run done.  Checkpoint: {os.path.basename(ckpt_file)}")

# ── Step B: restart from checkpoint → tend ────────────────────────────────────
print("\n" + "=" * 60)
print("Step B: load checkpoint → restart from day 1 to tend")
print("=" * 60)

s_restart = make_solver()
s_restart.set_checkpoint_prefix(ckpt_prefix)
s_restart.load_checkpoint(ckpt_file)
print(f"Restart time: {s_restart.restart_time()/86400:.2f} days")
s_restart.run(tend, dtout)

t_restart   = np.array(s_restart.time_points())
T_restart   = np.array(s_restart.peak_fuel_temperature())
gp_restart  = np.array(s_restart.gas_pressure())
bup_restart = np.array(s_restart.burnup())
fgr_restart = np.array(s_restart.fg_released())

print(f"\nRestart run: {len(t_restart)} points, "
      f"t=[{t_restart[0]/86400:.2f}, {t_restart[-1]/86400:.2f}] days")

# ── Verify against reference ──────────────────────────────────────────────────
ref_file = os.path.join(WORK_DIR, "ref_times.npy")
if not os.path.exists(ref_file):
    print("\n[INFO] Reference data not found — run script 0 first for comparison.")
else:
    t_ref   = np.load(os.path.join(WORK_DIR, "ref_times.npy"))
    T_ref   = np.load(os.path.join(WORK_DIR, "ref_T.npy"))
    gp_ref  = np.load(os.path.join(WORK_DIR, "ref_gp.npy"))
    bup_ref = np.load(os.path.join(WORK_DIR, "ref_bup.npy"))
    fgr_ref = np.load(os.path.join(WORK_DIR, "ref_fgr.npy"))

    print("\n" + "=" * 60)
    print("Verification: checkpoint restart vs reference (script 0)")
    print("=" * 60)

    tol_t = 1.0
    matches = [(i, j, t) for i, t in enumerate(t_ref)
               for j, tr in enumerate(t_restart) if abs(t - tr) < tol_t]
    print(f"Overlapping time points: {len(matches)}")

    def check(name, a_ref, a_rst, rtol=1e-6, atol=0.0):
        errs = [abs(float(a_ref[i]) - float(a_rst[j])) /
                max(abs(float(a_ref[i])), abs(float(a_rst[j])), atol, 1e-30)
                for i, j, _ in matches]
        mx = max(errs) if errs else 0.0
        ok = mx < rtol
        print(f"  {name:<35} max_rel_err={mx:.2e}  [{'PASS' if ok else 'FAIL'}]")
        return ok

    ok = True
    ok &= check("peak_fuel_temperature [K]", T_ref,   T_restart)
    ok &= check("gas_pressure [MPa]",        gp_ref,  gp_restart, rtol=1e-5)
    ok &= check("burnup [MWd/kgU]",          bup_ref, bup_restart)
    ok &= check("fg_released [mol]",         fgr_ref, fgr_restart, atol=1e-10)

    print()
    print("OVERALL:", "PASS" if ok else "FAIL",
          "— checkpoint restart reproduces the reference exactly." if ok
          else "— mismatch detected.")
