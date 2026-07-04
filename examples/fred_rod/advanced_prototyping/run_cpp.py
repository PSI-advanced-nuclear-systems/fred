"""
FRED-ROD advanced prototyping — Step 2: compiled C++ equivalent.

Uses fred.HeKrXe — the compiled C++ implementation of the same burnup-dependent
gap conductance model that was prototyped in run_python.py.

The C++ class (src/apps/fred_rod/gapmaterial/HeKrXe.hpp) implements identical
physics: Waltar-Reynolds FGR fraction, molar balance, geometric-mean conductivity.
After compilation it is available as fred.HeKrXe(bup_atpct=...) with no Python
overhead in the inner loop.

This file:
  1. Verifies that fred.HeKrXe conductance values match the Python prototype.
  2. Runs the same three simulations as run_python.py to confirm numerical
     agreement and demonstrate the cleaner API.

See README.md for compilation instructions.
"""

import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'build'))

import fred_rod as fred

# ===========================================================================
# Verify conductance agreement with the Python prototype
# ===========================================================================

# Inline reference matching HeKrXePy for spot-check (no import dependency)
def _py_k_gap(bup_atpct, T):
    bup_MWd = bup_atpct * 9.4
    fgr = 0.0
    if bup_MWd > 0.0:
        a = 4.7 / bup_MWd * (1.0 - math.exp(-bup_MWd / 5.9))
        fgr = max(0.0, 1.0 - a)
    n_He0 = 0.1e6 * 1.5e-6 / (8.314 * 293.15)
    fggen = 0.25 * (bup_atpct / 100.0) * 0.010 / 0.238
    fgrel = fgr * fggen
    n_Xe   = 0.8846 * fgrel; n_Kr = 0.0769 * fgrel; n_He_FG = 0.0385 * fgrel
    n_tot  = n_He0 + n_He_FG + n_Kr + n_Xe
    y_He = (n_He0 + n_He_FG) / n_tot; y_Kr = n_Kr / n_tot; y_Xe = n_Xe / n_tot
    k_He = 2.639e-3 * T**0.7085; k_Kr = 8.247e-5 * T**0.8363; k_Xe = 4.351e-5 * T**0.8618
    return k_He**y_He * k_Kr**y_Kr * k_Xe**y_Xe

print("=" * 64)
print("Python vs C++ conductance check")
print("=" * 64)
print(f"  {'bup [at%]':>10}  {'T [K]':>6}  {'Python k':>12}  {'C++ k':>12}  {'rel err':>10}")
print("  " + "-" * 55)

all_ok = True
for bup in [0.0, 1.0, 5.0]:
    for T in [500.0, 900.0, 1400.0]:
        k_py  = _py_k_gap(bup, T)
        k_cpp = fred.HeKrXe(bup).gap_conductivity(T)
        err   = abs(k_py - k_cpp) / max(abs(k_py), 1e-30)
        ok    = err < 1.0e-12
        all_ok &= ok
        print(f"  {bup:>10.1f}  {T:>6.0f}  {k_py:>12.6f}  {k_cpp:>12.6f}"
              f"  {err:>10.2e}  {'OK' if ok else 'FAIL'}")

print(f"\n  Agreement check: {'PASS' if all_ok else 'FAIL'}")

# ===========================================================================
# Geometry — identical to run_python.py
# ===========================================================================

def make_geometry():
    geom = fred.FuelRodGeometry()
    geom.nf   = 4
    geom.nc   = 3
    geom.nz   = 1
    geom.rfi0 = 0.0
    geom.rfo0 = 4.2e-3
    geom.rci0 = 4.3e-3
    geom.rco0 = 4.9e-3
    geom.dz0  = [0.10]
    geom.vgp  = 1.5e-6
    geom.ruff = 5.0e-6
    geom.rufc = 5.0e-6
    geom.build()
    return geom

T_COOLANT = 620.0
POWER     = 3.0e8

def make_solver(geom, gap_mat):
    s = fred.FredRodSolver(geom, fred.UO2(), fred.AIM1(), gap_mat)
    s.set_enable_heat_conduction(True)
    s.set_enable_stress_strain(True)
    s.set_initial_temperature(T_COOLANT)
    s.set_coolant_temperature([0.0, 2000.0], [T_COOLANT, T_COOLANT])
    s.set_power_density_history([0.0, 50.0, 2000.0], [0.0, POWER, POWER])
    s.set_coolant_pressure(0.5)
    s.set_initial_gas_pressure(0.1)
    s.set_tolerances(1e-6, 1e-8)
    return s

# ===========================================================================
# Run with fred.HeKrXe at three burnup levels
# ===========================================================================

CASES = [
    (0.0, "as-fabricated (pure He)"),
    (1.0, "moderate burnup (1 at%)"),
    (5.0, "high burnup    (5 at%)"),
]

geom  = make_geometry()
tend  = 1500.0
dtout = 150.0

print()
print("=" * 72)
print(f"Steady-state results with fred.HeKrXe  (Q = {POWER/1e6:.0f} MW/m³, "
      f"T_cool = {T_COOLANT:.0f} K)")
print("=" * 72)
print(f"  {'Burnup':>22}  {'T_peak [K]':>11}  {'gap [µm]':>10}  "
      f"{'sigh_clad [MPa]':>16}")
print("  " + "-" * 66)

for bup, label in CASES:
    gap = fred.HeKrXe(bup)
    s   = make_solver(geom, gap)
    s.run(tend, dtout)
    T_pk = s.peak_fuel_temperature()[-1]
    gw   = s.gap_width()[-1]
    sigh = s.clad_outer_hoop_stress()[-1]
    print(f"  {label:>22}  {T_pk:>11.1f}  {gw*1e6:>10.2f}  {sigh:>16.3f}")

print("""
Notes:
  Results should match run_python.py to within IDA solver tolerance (~1e-6).
  Any discrepancy indicates a divergence between the Python and C++ model
  implementations and should be investigated.

  fred.HeKrXe is a frozen snapshot of the mole-fraction calculation.
""")
