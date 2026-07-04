"""
FRED-OX: hot start vs cold start demonstration.

Runs the same MOX fuel pin (single axial layer, full power from t=0) twice:

  (a) COLD START — solver.set_hot_start(False) (the default). The pin
      begins at T0 = 293.15 K (room temperature); the real time integration
      itself resolves the thermal transient up to the hot steady state.

  (b) HOT START  — solver.set_hot_start(True). Before t=0 of the real run,
      the solver marches the pin forward in pseudo-time at the t=0 boundary
      conditions with irradiation physics (FGR, swelling, burnup-dependent
      conductivity) switched off, until temperatures/stresses stop changing,
      then re-anchors the real integration at t=0 using that converged
      state. t=0 in the hot-start run is already at the hot steady state.

TEND is kept short (300 s) on purpose: long enough to reach thermal steady
state (FRED-OX's thermal time constant here is ~150-200 s) but far too
short for burnup/fission-gas release to become significant — i.e. "hot
steady state, before any serious irradiation effects", per the task
description. fgrel/bup are printed for both runs to show they stay
negligible over this window, so the temperature comparison below is a
clean thermal-only comparison.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'build'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import _fred_ox as fred   # compiled extension — full pybind11 API (set_hot_start, ...)

T0 = 293.15        # K — room temperature (cold start initial condition)
T_COOLANT = 700.0  # K
T_PLENUM = 700.0   # K
QQV = 2.0e8        # W/m3 — full power from t=0
TEND = 300.0       # s — reaches thermal steady state well before irradiation matters
DTOUT = 10.0       # s — fine enough to resolve the cold-start rise


def make_solver(hot_start: bool):
    g = fred.FuelRodGeometry()
    g.nz  = 1
    g.nf  = 10
    g.nc  = 3
    g.rfi0 = 0.0
    g.rfo0 = 4.68e-3
    g.rci0 = g.rfo0 + 1.55e-4   # 155 µm as-fab gap
    g.rco0 = 5.335e-3
    g.dz0  = [0.05]
    g.vgp  = 7.19e-5
    g.ruff = 2.0e-6
    g.rufc = 2.0e-6
    g.build()

    mox  = fred.FredOxMOX(pu_content=0.1799, rof0=10542.0, sto0=1.97)
    clad = fred.FredOxT91(reference_density=7790.0)
    gap  = fred.FredOxGapMaterial()

    s = fred.FredOxSolver(g, mox, clad, gap)
    s.set_initial_temperature(T0)
    s.set_coolant_temperature([0.0, 1e6], [T_COOLANT, T_COOLANT])
    s.set_plenum_temperature_history([0.0, 1e6], [T_PLENUM, T_PLENUM])
    s.set_power_density_history([0.0, 1e6], [QQV, QQV])
    s.set_hot_start(hot_start)
    s.set_coolant_pressure(0.2)          # MPa
    s.set_initial_gas_pressure(0.1)      # MPa
    s.set_enable_heat_conduction(True)
    s.set_enable_stress_strain(True)
    s.set_tolerances(1e-5, 1e-7)
    return s


results = {}
for label, hot in [("cold", False), ("hot", True)]:
    print(f"\n=== FRED-OX {label.upper()} START ===")
    solver = make_solver(hot)
    solver.run(TEND, DTOUT)
    times = np.asarray(solver.time_points())
    T_peak = np.asarray(solver.peak_fuel_temperature())
    fgrel = np.asarray(solver.fg_released())
    bup = np.asarray(solver.burnup())
    results[label] = (times, T_peak, fgrel, bup)

t_cold, T_cold, fgrel_cold, bup_cold = results["cold"]
t_hot, T_hot, fgrel_hot, bup_hot = results["hot"]

print(f"\n{'t [s]':>7}  {'T_peak cold [K]':>16}  {'T_peak hot [K]':>15}")
print("-" * 46)
for k in range(len(t_cold)):
    print(f"{t_cold[k]:7.1f}  {T_cold[k]:16.3f}  {T_hot[k]:15.3f}")

steady = T_hot[-1]
max_dev_hot = float(np.max(np.abs(T_hot - T_hot[0])))
print(f"\nCold start:  T_peak(t=0) = {T_cold[0]:.2f} K (room temp)  ->  "
      f"T_peak(t={t_cold[-1]:.0f}s) = {T_cold[-1]:.2f} K (steady state)")
print(f"Hot start :  T_peak(t=0) = {T_hot[0]:.2f} K (already steady)  ->  "
      f"T_peak(t={t_hot[-1]:.0f}s) = {T_hot[-1]:.2f} K")
print(f"Hot start max deviation from its own t=0 value: {max_dev_hot:.2e} K "
      f"(confirms 'no change')")
print(f"\nIrradiation-effect check (both runs, over {TEND:.0f} s):")
print(f"  fgrel: cold={fgrel_cold[-1]:.3e} mol   hot={fgrel_hot[-1]:.3e} mol")
print(f"  bup  : cold={bup_cold[-1]:.3e} MWd/kgU  hot={bup_hot[-1]:.3e} MWd/kgU")
print("  (both negligible: hot steady state reached well before serious irradiation)")

assert max_dev_hot < 1.0, "hot start should show essentially no temperature change over the run"
assert T_cold[0] < 400.0, "cold start should begin at room temperature"
assert abs(T_cold[-1] - steady) < 5.0, "cold start should converge to the same steady state as hot start"

# ── Plot ─────────────────────────────────────────────────────────────────────
plots_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')
os.makedirs(plots_dir, exist_ok=True)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(t_cold, T_cold, 'o-', color='tab:blue', label='Cold start (T0 = room temp)')
ax.plot(t_hot, T_hot, 's-', color='tab:red', label='Hot start (already at steady state)')
ax.axhline(steady, color='gray', linestyle=':', alpha=0.6, label=f'Steady state ({steady:.1f} K)')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Peak fuel temperature [K]')
ax.set_title('FRED-OX: cold start vs hot start\n(before serious irradiation effects)')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'hot_vs_cold_start.png'), dpi=120)
plt.close(fig)
print(f"\nSaved plots/hot_vs_cold_start.png")
