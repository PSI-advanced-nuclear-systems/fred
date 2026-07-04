"""
FRED-OX restart_snapshot — Script 2: Continuation via snapshot.

Two-phase irradiation / transient workflow:

  Phase 1 (base irradiation): run from t0 to t_half.
    → Snapshot saved automatically at t_half (ox_frame1.snapshot).
    → This represents the fully irradiated fuel state after day 1.

  Phase 2 (transient / extended irradiation): load the snapshot.
    → Simulation time resets to 0.0.
    → Physical state (temperatures, burnup, fission gas, strains) is preserved.
    → Run from 0 to (tend - t_half).

Combined trajectory (offset phase-2 times by t_half) should match the
reference full run (script 0).
"""

import sys, os
import numpy as np

BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import fred_ox as fred

WORK_DIR = os.path.join(os.path.dirname(__file__), 'ckpt')
os.makedirs(WORK_DIR, exist_ok=True)

tend     = 3.0 * 86400
dtout    = 86400.0
t_half   = 86400.0
t_phase2 = tend - t_half   # duration of phase 2

snap_prefix = os.path.join(WORK_DIR, "ox")
snap_file   = os.path.join(WORK_DIR, "ox_frame1.snapshot")


def make_solver(t_end_bc=None):
    if t_end_bc is None:
        t_end_bc = tend
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
    t_bc = [0.0, max(t_end_bc, tend)]
    s.set_power_density_history(t_bc, [1.5e8, 1.5e8])
    s.set_coolant_temperature(  t_bc, [620.0, 620.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(620.0)
    s.set_initial_gas_pressure(0.1)
    s.set_tolerances(1.0e-5, 1.0e-7)
    return s


# ── Phase 1: run 0 → t_half ──────────────────────────────────────────────────
print("=" * 60)
print(f"Phase 1: irradiation 0 → day 1  (snapshot saved at end)")
print("=" * 60)

s1 = make_solver()
s1.set_snapshot_prefix(snap_prefix)   # final snapshot always saved
s1.run(t_half, dtout)

t1   = np.array(s1.time_points())
T1   = np.array(s1.peak_fuel_temperature())
gp1  = np.array(s1.gas_pressure())
bup1 = np.array(s1.burnup())
fgr1 = np.array(s1.fg_released())

print(f"Phase 1 done: {len(t1)} points, t=[{t1[0]/86400:.2f}, {t1[-1]/86400:.2f}] d")
print(f"  Burnup at day 1:      {bup1[-1]:.4f} MWd/kgU")
print(f"  Snapshot: {os.path.basename(snap_file)}")

# ── Phase 2: load snapshot, run 0 → t_phase2 ─────────────────────────────────
print("\n" + "=" * 60)
print(f"Phase 2: load snapshot → run 0 → {t_phase2/86400:.1f} days (time reset)")
print("=" * 60)

s2 = make_solver(t_end_bc=t_phase2)
s2.load_snapshot(snap_file)
print(f"  restart_time = {s2.restart_time():.1f} s  (0 = time reset confirmed)")
print(f"  (burnup/fission-gas state preserved from phase 1)")

s2.run(t_phase2, dtout)

t2   = np.array(s2.time_points())
T2   = np.array(s2.peak_fuel_temperature())
gp2  = np.array(s2.gas_pressure())
bup2 = np.array(s2.burnup())
fgr2 = np.array(s2.fg_released())

print(f"Phase 2 done: {len(t2)} points, t=[{t2[0]/86400:.2f}, {t2[-1]/86400:.2f}] d")

# ── Combined ──────────────────────────────────────────────────────────────────
t_combined   = np.concatenate([t1, t2[1:] + t_half])
T_combined   = np.concatenate([T1, T2[1:]])
gp_combined  = np.concatenate([gp1, gp2[1:]])
bup_combined = np.concatenate([bup1, bup2[1:]])
fgr_combined = np.concatenate([fgr1, fgr2[1:]])

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

    print()
    print("OVERALL:", "PASS" if ok else "FAIL",
          "— continuation matches reference exactly." if ok
          else "— mismatch detected.")
