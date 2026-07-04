"""
FRED-OX: Checkpoint save/load restart example

Demonstrates:
  1. Reference full run (0 → tend) with checkpoints.
  2. First-half run (0 → t_half), then restart from checkpoint to tend.
  3. Verify that the restarted run matches the reference at overlapping times.

Checkpoint file: <prefix>_checkpoint.chk  (single file, overwritten each step)

Geometry: MOX pin / T91 cladding / He-like gap, single axial layer.
"""

import sys, os
import numpy as np

BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import fred_ox as fred

WORK_DIR = os.path.join(os.path.dirname(__file__), 'ckpt')
os.makedirs(WORK_DIR, exist_ok=True)


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

    tend = 3.0 * 86400

    s = fred.FredOxSolver(g, mox, clad, gap)
    s.set_power_density_history([0.0, tend], [1.5e8, 1.5e8])
    s.set_coolant_temperature(  [0.0, tend], [620.0, 620.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(620.0)
    s.set_initial_gas_pressure(0.1)
    s.set_tolerances(1.0e-5, 1.0e-7)
    return s


tend   = 3.0 * 86400   # 3 days
dtout  = 86400.0        # output daily
t_half = 86400.0        # stop first run here

ckpt_prefix = os.path.join(WORK_DIR, "ox")
ckpt_file   = os.path.join(WORK_DIR, "ox_checkpoint.chk")

# ── Reference: full run ───────────────────────────────────────────────────────
print("=" * 60)
print("Reference: full run 0 → {:.0f} days".format(tend / 86400))
print("=" * 60)

s_ref = make_solver()
s_ref.run(tend, dtout)

t_ref  = np.array(s_ref.time_points())
T_ref  = np.array(s_ref.peak_fuel_temperature())
gp_ref = np.array(s_ref.gas_pressure())
bup_ref = np.array(s_ref.burnup())
fgr_ref = np.array(s_ref.fg_released())

print(f"Full run: {len(t_ref)} output points, t=[{t_ref[0]/86400:.2f}, {t_ref[-1]/86400:.2f}] d")

# ── Case 1: first-half run (checkpoint saved at t_half) ───────────────────────
print("\n" + "=" * 60)
print("Case 1: first-half run 0 → day 1  (checkpoints enabled)")
print("=" * 60)

s1 = make_solver()
s1.set_checkpoint_prefix(ckpt_prefix)
s1.run(t_half, dtout)

print(f"First-half done.  Checkpoint: {os.path.basename(ckpt_file)}")

# ── Case 2: restart from checkpoint at t_half, continue to tend ───────────────
print("\n" + "=" * 60)
print("Case 2: restart from checkpoint at day 1 → tend")
print("=" * 60)

s2 = make_solver()
s2.set_checkpoint_prefix(ckpt_prefix)
s2.load_checkpoint(ckpt_file)
print(f"Restart time: {s2.restart_time()/86400:.2f} days")
s2.run(tend, dtout)

t_restart   = np.array(s2.time_points())
T_restart   = np.array(s2.peak_fuel_temperature())
gp_restart  = np.array(s2.gas_pressure())
bup_restart = np.array(s2.burnup())
fgr_restart = np.array(s2.fg_released())

print(f"\nRestart run: {len(t_restart)} output points, "
      f"t=[{t_restart[0]/86400:.2f}, {t_restart[-1]/86400:.2f}] d")

# ── Case 3: verify restart matches reference ──────────────────────────────────
print("\n" + "=" * 60)
print("Case 3: verify restart == reference at overlapping time points")
print("=" * 60)

tol_t = 1.0  # s
matches = [(i, j, t) for i, t in enumerate(t_ref)
           for j, tr in enumerate(t_restart) if abs(t - tr) < tol_t]
print(f"Overlapping output points: {len(matches)}")

def check_match(name, arr_ref, arr_restart, rtol=1.0e-6, atol=0.0):
    errs = []
    for i_f, i_r, _ in matches:
        v_f = float(arr_ref[i_f]);  v_r = float(arr_restart[i_r])
        ref = max(abs(v_f), abs(v_r), atol, 1.0e-30)
        errs.append(abs(v_f - v_r) / ref)
    mx = max(errs) if errs else 0.0
    ok = mx < rtol
    print(f"  {name:<35} max_rel_err={mx:.2e}  [{'PASS' if ok else 'FAIL'}]")
    return ok

all_pass = True
all_pass &= check_match("peak_fuel_temperature [K]", T_ref,   T_restart)
all_pass &= check_match("gas_pressure [MPa]",        gp_ref,  gp_restart, rtol=1.0e-5)
all_pass &= check_match("burnup [MWd/kgU]",          bup_ref, bup_restart)
all_pass &= check_match("fg_released [mol]",         fgr_ref, fgr_restart, atol=1.0e-10)

print()
print("OVERALL:", "PASS" if all_pass else "FAIL",
      "— checkpoint restart reproduces full run exactly." if all_pass
      else "— mismatch detected.")
