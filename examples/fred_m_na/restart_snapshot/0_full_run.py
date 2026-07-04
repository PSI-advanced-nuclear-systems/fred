"""
FRED-M-Na restart_snapshot — Script 0: Full simulation from t0 to tend.

Runs the complete irradiation with:
  • Checkpointing enabled (mna_checkpoint.chk overwritten each output step).
  • A snapshot at t_half (user-specified) and at tend (always automatic).

This is the reference run.  The checkpoint and snapshots it produces are
used by scripts 1 and 2 respectively.

Geometry: U-19Pu-10Zr / HT-9, EBR-II-like pin, 3 days at steady power.
"""

import sys, os
import numpy as np

BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import fred_m_na as fred

WORK_DIR = os.path.join(os.path.dirname(__file__), 'ckpt')
os.makedirs(WORK_DIR, exist_ok=True)

# ── Simulation parameters ─────────────────────────────────────────────────────
tend   = 3.0 * 86400   # 3 days [s]
dtout  = 86400.0        # output every day [s]
t_half = 86400.0        # mid-point snapshot time (day 1)


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
    # Slow ramp avoids IDACalcIC failure; full power from t=100 s onward.
    s.set_power_density_history([0.0, 50.0, 100.0, tend], [0.0, 0.0, qqv, qqv])
    s.set_coolant_temperature([0.0, tend], [643.0, 643.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(673.0)
    s.set_initial_gas_pressure(0.1)
    return s


# ── Full run ──────────────────────────────────────────────────────────────────
print("=" * 60)
print("FRED-M-Na: full run  t=0 → {:.0f} days".format(tend / 86400))
print("=" * 60)

s = make_solver()

# Checkpoint: fault-recovery file (overwritten each output step).
s.set_checkpoint_prefix(os.path.join(WORK_DIR, "mna"))

# Snapshot at t_half; the final snapshot (at tend) is always saved automatically.
s.set_snapshot_prefix(os.path.join(WORK_DIR, "mna"),
                      snapshot_times=[t_half])

s.run(tend, dtout)

t_out   = np.array(s.times())
T_peak  = np.array(s.peak_fuel_temperature())
gp_out  = np.array(s.gas_pressure())
bup_out = np.array(s.burnup())
fgr_out = np.array(s.fg_released())
sw_out  = np.array(s.grsis_swelling_total())
gap_out = np.array(s.gap_width())

print(f"\nOutput: {len(t_out)} time points, "
      f"t=[{t_out[0]/86400:.2f}, {t_out[-1]/86400:.2f}] days")
print(f"  Peak fuel temperature (final): {T_peak[-1]:.2f} K")
print(f"  Gas pressure (final):          {gp_out[-1]:.4f} MPa")
print(f"  Burnup (final):                {bup_out[-1]:.4f} MWd/kgU")
print(f"  GRSIS swelling (final):        {sw_out[-1]:.6f}")

# Save reference results for comparison in scripts 1 and 2.
np.save(os.path.join(WORK_DIR, "ref_times.npy"),  t_out)
np.save(os.path.join(WORK_DIR, "ref_T.npy"),      T_peak)
np.save(os.path.join(WORK_DIR, "ref_gp.npy"),     gp_out)
np.save(os.path.join(WORK_DIR, "ref_bup.npy"),    bup_out)
np.save(os.path.join(WORK_DIR, "ref_fgr.npy"),    fgr_out)
np.save(os.path.join(WORK_DIR, "ref_sw.npy"),     sw_out)
np.save(os.path.join(WORK_DIR, "ref_gap.npy"),    gap_out)

print("\nFiles written to", WORK_DIR)
print("  mna_checkpoint.chk      ← fault-recovery checkpoint (at tend)")
print("  mna_frame1.snapshot     ← snapshot at t_half (day 1)")
print("  mna_frame2.snapshot     ← snapshot at tend   (day 3, final state)")
print("  ref_*.npy               ← reference results for comparison")
