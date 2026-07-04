"""
FRED-ROD: hot start vs cold start demonstration.

Runs the *same* rod (dummy fuel/cladding/gap, full power from t=0) twice:

  (a) COLD START  — solver.set_hot_start(False) (the default). The rod
      begins at T0 = 293.15 K (room temperature) and the real time
      integration itself carries it through the transient up to the hot
      steady state. This is what "hot start" replaces: an actual multi-step
      thermal transient has to be resolved before the rod is at operating
      conditions.

  (b) HOT START   — solver.set_hot_start(True). Before t=0 of the real run,
      the solver marches the rod forward in pseudo-time at the t=0 boundary
      conditions (irradiation physics off — FRED-ROD has none anyway) until
      temperatures/stresses stop changing, then re-anchors the real
      integration at t=0 using that converged state. So even the very first
      output point (t=0) is already at the hot steady state, and every
      later point is (up to solver tolerance) identical to it — "no change".

No irradiation physics exists in FRED-ROD, so there is no "before serious
irradiation effects" caveat here (see the FRED-OX / FRED-M-Na hot_start_demo
examples for that side of the comparison) — this example is a pure
heat-conduction + stress-strain steady-state comparison.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'build'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import _fred_rod as fred   # compiled extension — full pybind11 API (set_hot_start, ...)

T0 = 293.15       # K — room temperature (cold start initial condition)
T_COOLANT = 700.0 # K
QQV = 1.0e8       # W/m3 — full power from t=0
TEND = 200.0      # s — well past the thermal time constant
DTOUT = 5.0       # s — fine enough to resolve the cold-start rise


def make_solver(hot_start: bool):
    g = fred.FuelRodGeometry()
    g.nf   = 3
    g.nc   = 2
    g.nz   = 1
    g.rfi0 = 0.0        # solid pellet
    g.rfo0 = 4.2e-3     # m
    g.rci0 = 4.3e-3     # m   (100 µm as-fab gap)
    g.rco0 = 4.9e-3     # m
    g.dz0  = [0.10]     # m
    g.vgp  = 1.0e-6     # m3
    g.ruff = 5.0e-6     # m
    g.rufc = 5.0e-6     # m
    g.build()

    fuel = fred.DummyFuelPellet()
    clad = fred.DummyCladding()
    gap  = fred.DummyGapMaterial()

    s = fred.FredRodSolver(g, fuel, clad, gap)
    s.set_enable_heat_conduction(True)
    s.set_enable_stress_strain(True)
    s.set_initial_temperature(T0)
    s.set_coolant_temperature([0.0, 1e6], [T_COOLANT, T_COOLANT])
    s.set_power_density_history([0.0, 1e6], [QQV, QQV])
    s.set_hot_start(hot_start)
    s.set_coolant_pressure(5.0)          # MPa
    s.set_initial_gas_pressure(0.1)      # MPa
    s.set_tolerances(1e-6, 1e-8)
    return s


results = {}
for label, hot in [("cold", False), ("hot", True)]:
    print(f"\n=== FRED-ROD {label.upper()} START ===")
    solver = make_solver(hot)
    solver.run(TEND, DTOUT)
    times = np.asarray(solver.time_points())
    T = np.asarray(solver.temperatures())   # [n_steps, nf+nc]
    results[label] = (times, T)

t_cold, T_cold = results["cold"]
t_hot, T_hot = results["hot"]

nf, nc = 3, 2
T_ctr_cold = T_cold[:, 0]     # fuel centreline
T_ctr_hot = T_hot[:, 0]

print(f"\n{'t [s]':>7}  {'T_ctr cold [K]':>15}  {'T_ctr hot [K]':>15}")
print("-" * 45)
for k in range(len(t_cold)):
    print(f"{t_cold[k]:7.1f}  {T_ctr_cold[k]:15.3f}  {T_ctr_hot[k]:15.3f}")

steady = T_ctr_hot[-1]
max_dev_hot = float(np.max(np.abs(T_ctr_hot - T_ctr_hot[0])))
print(f"\nCold start:  T_ctr(t=0) = {T_ctr_cold[0]:.2f} K (room temp)  ->  "
      f"T_ctr(t={t_cold[-1]:.0f}s) = {T_ctr_cold[-1]:.2f} K (steady state)")
print(f"Hot start :  T_ctr(t=0) = {T_ctr_hot[0]:.2f} K (already steady)  ->  "
      f"T_ctr(t={t_hot[-1]:.0f}s) = {T_ctr_hot[-1]:.2f} K")
print(f"Hot start max deviation from its own t=0 value: {max_dev_hot:.2e} K "
      f"(confirms 'no change')")
assert max_dev_hot < 0.1, "hot start should show essentially no temperature change over the run"
assert T_ctr_cold[0] < 400.0, "cold start should begin at room temperature"
assert abs(T_ctr_cold[-1] - steady) < 1.0, "cold start should converge to the same steady state as hot start"

# ── Plot ─────────────────────────────────────────────────────────────────────
plots_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')
os.makedirs(plots_dir, exist_ok=True)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(t_cold, T_ctr_cold, 'o-', color='tab:blue', label='Cold start (T0 = room temp)')
ax.plot(t_hot, T_ctr_hot, 's-', color='tab:red', label='Hot start (already at steady state)')
ax.axhline(steady, color='gray', linestyle=':', alpha=0.6, label=f'Steady state ({steady:.1f} K)')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Fuel centreline temperature [K]')
ax.set_title('FRED-ROD: cold start vs hot start')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'hot_vs_cold_start.png'), dpi=120)
plt.close(fig)
print(f"\nSaved plots/hot_vs_cold_start.png")
