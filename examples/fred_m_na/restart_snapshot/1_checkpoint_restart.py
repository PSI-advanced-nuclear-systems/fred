"""
FRED-M-Na restart_snapshot — Script 1: Restart from checkpoint.

Scenario: simulation crashed after day 1.  The checkpoint (mna_checkpoint.chk)
contains the full fuel state (temperatures, Zr distribution, GRSIS bubbles,
burnup, fission gas) at the end of day 1.  This script loads that checkpoint
and resumes the simulation to tend.

Simulation time is preserved (restarts from t_half, NOT reset to 0).

Run script 0 first to generate the reference data and checkpoint file,
or let Step A in this script generate the checkpoint automatically.
"""

import sys, os
import numpy as np

BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import fred_m_na as fred

WORK_DIR = os.path.join(os.path.dirname(__file__), 'ckpt')
os.makedirs(WORK_DIR, exist_ok=True)

# ── Simulation parameters (must match script 0) ───────────────────────────────
tend   = 3.0 * 86400
dtout  = 86400.0
t_half = 86400.0   # checkpoint restart point (day 1)

ckpt_prefix = os.path.join(WORK_DIR, "mna")
ckpt_file   = os.path.join(WORK_DIR, "mna_checkpoint.chk")


def make_solver():
    g = fred.FuelRodGeometry()
    g.nf  = 4;  g.nc  = 3;  g.nz  = 1
    g.rfi0 = 0.0;  g.rfo0 = 2.35e-3
    g.rci0 = 2.45e-3;  g.rco0 = 2.92e-3
    g.dz0  = [0.343]
    g.vgp  = 2.5e-6;  g.ruff = 5e-6;  g.rufc = 5e-6
    g.build()

    fuel = fred.UPuZr(pu_weight_frac=0.19, zr_weight_frac=0.10,
                      reference_density=15700.0)
    clad = fred.HT9(reference_density=7750.0)

    s = fred.FredMNaSolver(g, fuel, clad)
    qqv = 3.67e9
    s.set_power_density_history([0.0, 50.0, 100.0, tend], [0.0, 0.0, qqv, qqv])
    s.set_coolant_temperature([0.0, tend], [643.0, 643.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(673.0)
    s.set_initial_gas_pressure(0.1)
    return s


# ── Step A: run 0 → t_half (simulates the partial run before the crash) ───────
print("=" * 60)
print(f"Step A: partial run 0 → day 1  (checkpoints enabled)")
print("=" * 60)

s_partial = make_solver()
s_partial.set_checkpoint_prefix(ckpt_prefix)
s_partial.run(t_half, dtout)

print(f"Partial run done.  Checkpoint: {os.path.basename(ckpt_file)}")

# ── Step B: load checkpoint at t_half, continue to tend ───────────────────────
print("\n" + "=" * 60)
print("Step B: load checkpoint → restart from day 1 to tend")
print("=" * 60)

s_restart = make_solver()
s_restart.set_checkpoint_prefix(ckpt_prefix)   # continue writing checkpoints
s_restart.load_checkpoint(ckpt_file)
print(f"Loaded checkpoint.  Restart time: {s_restart.restart_time()/86400:.2f} days")

s_restart.run(tend, dtout)

t_restart   = np.array(s_restart.times())
T_restart   = np.array(s_restart.peak_fuel_temperature())
gp_restart  = np.array(s_restart.gas_pressure())
bup_restart = np.array(s_restart.burnup())
fgr_restart = np.array(s_restart.fg_released())
sw_restart  = np.array(s_restart.grsis_swelling_total())

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
    sw_ref  = np.load(os.path.join(WORK_DIR, "ref_sw.npy"))

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
    ok &= check("grsis_swelling_total [-]",  sw_ref,  sw_restart,  atol=1e-12)

    print()
    print("OVERALL:", "PASS" if ok else "FAIL",
          "— checkpoint restart reproduces the reference exactly." if ok
          else "— mismatch detected.")
