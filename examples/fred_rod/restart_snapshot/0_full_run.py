"""
FRED-ROD restart_snapshot — Script 0: Full simulation from t0 to tend.

Runs the complete simulation with:
  • Checkpointing enabled (checkpoint.chk overwritten each output step).
  • A snapshot saved at t_half (user-specified) and always at tend (automatic).

This script is the reference run.  The checkpoint and snapshots it produces
are loaded by scripts 1 and 2 respectively.

Geometry: UO2 pin / AIM1 cladding / He gap, 2 axial layers.
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
tend   = 5000.0   # total simulation time [s]
dtout  = 1000.0   # output interval [s]
t_half = 3000.0   # mid-point for snapshot and checkpoint restart


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


# ── Full run ──────────────────────────────────────────────────────────────────
print("=" * 60)
print("FRED-ROD: full run  t=0 → {:.0f} s".format(tend))
print("=" * 60)

s = make_solver()

# Checkpoint: fault-recovery file (overwritten each step).
s.set_checkpoint_prefix(os.path.join(WORK_DIR, "rod"))

# Snapshot at t_half; final snapshot (at tend) is always saved automatically.
s.set_snapshot_prefix(os.path.join(WORK_DIR, "rod"),
                      snapshot_times=[t_half])

s.run(tend, dtout)

t_out   = np.array(s.time_points())
T_peak  = np.array(s.peak_fuel_temperature())
gap_out = np.array(s.gap_width())
rco_out = np.array(s.clad_outer_radius())

print(f"\nOutput: {len(t_out)} time points, "
      f"t=[{t_out[0]:.0f}, {t_out[-1]:.0f}] s")
print(f"  Peak fuel temperature (final): {T_peak[-1]:.2f} K")
print(f"  Gap width (final):             {gap_out[-1]*1e6:.2f} µm")

# Save results for comparison in the other scripts.
np.save(os.path.join(WORK_DIR, "ref_times.npy"),  t_out)
np.save(os.path.join(WORK_DIR, "ref_T.npy"),      T_peak)
np.save(os.path.join(WORK_DIR, "ref_gap.npy"),    gap_out)
np.save(os.path.join(WORK_DIR, "ref_rco.npy"),    rco_out)

print("\nFiles written to", WORK_DIR)
print("  rod_checkpoint.chk      ← checkpoint at tend (fault recovery)")
print("  rod_frame1.snapshot     ← snapshot at t_half (for continuation)")
print("  rod_frame2.snapshot     ← snapshot at tend   (final state)")
print("  ref_*.npy               ← reference results for comparison")
