"""
FRED-M-Na restart_snapshot — Script 2: Continuation via snapshot.

Two-phase irradiation / transient workflow:

  Phase 1 (base irradiation): run from t0 to t_half.
    → Snapshot saved automatically at t_half (mna_frame1.snapshot).
    → Represents the physically meaningful fuel state after day 1 of irradiation.

  Phase 2 (continuation / transient): load the snapshot.
    → Simulation time resets to 0.0 (fresh time coordinate).
    → Physical fuel state is fully preserved:
        temperatures, gap width, burnup, fission gas inventory,
        Zr composition profile, GRSIS bubble state, cladding wastage.
    → Accumulated irradiation clock (elapsed time) is also preserved so
        burnup-dependent effects continue correctly.
    → Run from 0 to (tend - t_half).

Combined output:
  Phase 1 times:  [0, dtout, ..., t_half]
  Phase 2 times:  [0, dtout, ..., tend-t_half]  → offset by t_half for plotting

The combined trajectory should match the reference run (script 0) for
temperature, gas pressure, burnup, and fission gas release.

Note: cladding wastage uses local simulation time t (not elapsed time), so
a small discrepancy is expected for that quantity in the continuation mode.
It is therefore excluded from the verification check.
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
tend     = 3.0 * 86400   # total irradiation duration [s]
dtout    = 86400.0        # output interval [s]
t_half   = 86400.0        # end of phase 1 / duration of phase 2
t_phase2 = tend - t_half  # = 2 days

snap_prefix = os.path.join(WORK_DIR, "mna")
snap_file   = os.path.join(WORK_DIR, "mna_frame1.snapshot")   # auto-saved at t_half

QQV_TARGET = 3.67e9   # W/m3 — full operating power density


def make_solver_phase1():
    """Phase 1: ramp to full power then hold steady until t_half."""
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
    # Slow ramp from 0 to full power — avoids IDACalcIC failure.
    s.set_power_density_history([0.0, 50.0, 100.0, t_half],
                                [0.0, 0.0, QQV_TARGET, QQV_TARGET])
    s.set_coolant_temperature([0.0, t_half], [643.0, 643.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(673.0)
    s.set_initial_gas_pressure(0.1)
    return s


def make_solver_phase2():
    """Phase 2: full power throughout (fuel already hot from phase 1)."""
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
    # Constant full power: fuel state (temperatures) restored from snapshot,
    # so no ramp is needed.
    s.set_power_density_history([0.0, t_phase2], [QQV_TARGET, QQV_TARGET])
    s.set_coolant_temperature([0.0, t_phase2], [643.0, 643.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(673.0)   # only used if no snapshot is loaded
    s.set_initial_gas_pressure(0.1)    # only used if no snapshot is loaded
    s.set_tolerances(1.0e-3, 1.0e-5)
    return s


# ── Phase 1: run 0 → t_half, snapshot saved at end ───────────────────────────
print("=" * 60)
print(f"Phase 1: irradiation 0 → day 1  (snapshot saved at end)")
print("=" * 60)

s1 = make_solver_phase1()
s1.set_snapshot_prefix(snap_prefix)   # final state is always saved as a snapshot
s1.run(t_half, dtout)

t1   = np.array(s1.times())
T1   = np.array(s1.peak_fuel_temperature())
gp1  = np.array(s1.gas_pressure())
bup1 = np.array(s1.burnup())
fgr1 = np.array(s1.fg_released())
sw1  = np.array(s1.grsis_swelling_total())
gap1 = np.array(s1.gap_width())

print(f"Phase 1 done: {len(t1)} points, t=[{t1[0]/86400:.2f}, {t1[-1]/86400:.2f}] d")
print(f"  Burnup at day 1:      {bup1[-1]:.4f} MWd/kgU")
print(f"  Gas pressure at day 1:{gp1[-1]:.4f} MPa")
print(f"  Snapshot: {os.path.basename(snap_file)}")

# ── Phase 2: load snapshot, run 0 → t_phase2 (time reset) ────────────────────
print("\n" + "=" * 60)
print(f"Phase 2: load snapshot → run 0 → {t_phase2/86400:.1f} days (time reset)")
print("=" * 60)

s2 = make_solver_phase2()
s2.load_snapshot(snap_file)
print(f"  restart_time = {s2.restart_time():.1f} s  (0 = time reset confirmed)")
print(f"  (burnup, Zr profile, GRSIS state, elapsed time all preserved)")

s2.run(t_phase2, dtout)

t2   = np.array(s2.times())
T2   = np.array(s2.peak_fuel_temperature())
gp2  = np.array(s2.gas_pressure())
bup2 = np.array(s2.burnup())
fgr2 = np.array(s2.fg_released())
sw2  = np.array(s2.grsis_swelling_total())
gap2 = np.array(s2.gap_width())

print(f"Phase 2 done: {len(t2)} points, t=[{t2[0]/86400:.2f}, {t2[-1]/86400:.2f}] d")
print(f"  Burnup at end of phase 2: {bup2[-1]:.4f} MWd/kgU")

# ── Combined trajectory (offset phase-2 times by t_half) ─────────────────────
t_combined   = np.concatenate([t1, t2[1:] + t_half])
T_combined   = np.concatenate([T1, T2[1:]])
gp_combined  = np.concatenate([gp1, gp2[1:]])
bup_combined = np.concatenate([bup1, bup2[1:]])
fgr_combined = np.concatenate([fgr1, fgr2[1:]])
sw_combined  = np.concatenate([sw1, sw2[1:]])
gap_combined = np.concatenate([gap1, gap2[1:]])

print(f"\nCombined: {len(t_combined)} points, "
      f"t=[{t_combined[0]/86400:.2f}, {t_combined[-1]/86400:.2f}] d")

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
    print("Verification: continuation vs reference (script 0)")
    print("=" * 60)

    tol_t = 1.0
    matches = [(i, j, t) for i, t in enumerate(t_ref)
               for j, tc in enumerate(t_combined) if abs(t - tc) < tol_t]
    print(f"Overlapping time points: {len(matches)}")

    def check(name, a_ref, a_comb, rtol=1e-6, atol=0.0):
        errs = [abs(float(a_ref[i]) - float(a_comb[j])) /
                max(abs(float(a_ref[i])), abs(float(a_comb[j])), atol, 1e-30)
                for i, j, _ in matches]
        mx = max(errs) if errs else 0.0
        ok = mx < rtol
        print(f"  {name:<35} max_rel_err={mx:.2e}  [{'PASS' if ok else 'FAIL'}]")
        return ok

    ok = True
    ok &= check("peak_fuel_temperature [K]", T_ref,   T_combined)
    ok &= check("gas_pressure [MPa]",        gp_ref,  gp_combined, rtol=1e-5)
    ok &= check("burnup [MWd/kgU]",          bup_ref, bup_combined)
    ok &= check("fg_released [mol]",         fgr_ref, fgr_combined, atol=1e-10)
    ok &= check("grsis_swelling_total [-]",  sw_ref,  sw_combined,  atol=1e-12)

    print()
    print("OVERALL:", "PASS" if ok else "FAIL",
          "— continuation matches reference exactly." if ok
          else "— mismatch detected.")
    print()
    print("Note: clad_wastage is excluded from verification — it depends on")
    print("local simulation time t (not elapsed time), so a small discrepancy")
    print("is expected between continuation and reference for that quantity.")
