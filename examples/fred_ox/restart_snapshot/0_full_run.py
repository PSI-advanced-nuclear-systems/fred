"""
FRED-OX restart_snapshot — Script 0: Full simulation from t0 to tend.

Runs the complete irradiation with:
  • Checkpointing enabled (ox_checkpoint.chk overwritten each step).
  • A snapshot at t_half (user-specified) and at tend (always automatic).

This is the reference run.  Checkpoint and snapshots are used by scripts 1 and 2.

Geometry: MOX pin / T91 cladding / FredOxGapMaterial, 1 axial layer, 3 days.
"""

import sys, os
import numpy as np

BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import fred_ox as fred

WORK_DIR = os.path.join(os.path.dirname(__file__), 'ckpt')
os.makedirs(WORK_DIR, exist_ok=True)

# ── Simulation parameters ─────────────────────────────────────────────────────
tend   = 3.0 * 86400   # 3 days [s]
dtout  = 86400.0        # output every day [s]
t_half = 86400.0        # mid-point (day 1)


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


# ── Full run ──────────────────────────────────────────────────────────────────
print("=" * 60)
print("FRED-OX: full run  t=0 → {:.0f} days".format(tend / 86400))
print("=" * 60)

s = make_solver()
s.set_checkpoint_prefix(os.path.join(WORK_DIR, "ox"))
s.set_snapshot_prefix(os.path.join(WORK_DIR, "ox"),
                      snapshot_times=[t_half])

s.run(tend, dtout)

t_out  = np.array(s.time_points())
T_peak = np.array(s.peak_fuel_temperature())
gp_out = np.array(s.gas_pressure())
bup_out = np.array(s.burnup())
fgr_out = np.array(s.fg_released())

print(f"\nOutput: {len(t_out)} time points, "
      f"t=[{t_out[0]/86400:.2f}, {t_out[-1]/86400:.2f}] days")
print(f"  Peak fuel temperature (final): {T_peak[-1]:.2f} K")
print(f"  Gas pressure (final):          {gp_out[-1]:.4f} MPa")
print(f"  Burnup (final):                {bup_out[-1]:.4f} MWd/kgU")

np.save(os.path.join(WORK_DIR, "ref_times.npy"),  t_out)
np.save(os.path.join(WORK_DIR, "ref_T.npy"),      T_peak)
np.save(os.path.join(WORK_DIR, "ref_gp.npy"),     gp_out)
np.save(os.path.join(WORK_DIR, "ref_bup.npy"),    bup_out)
np.save(os.path.join(WORK_DIR, "ref_fgr.npy"),    fgr_out)

print("\nFiles in", WORK_DIR)
print("  ox_checkpoint.chk      ← fault-recovery checkpoint (at tend)")
print("  ox_frame1.snapshot     ← snapshot at t_half (day 1)")
print("  ox_frame2.snapshot     ← snapshot at tend   (day 3, final state)")
print("  ref_*.npy              ← reference results for comparison")
