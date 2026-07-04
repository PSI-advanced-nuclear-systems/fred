"""
FRED-M-Na: hot start vs cold start demonstration.

Runs the same U-19Pu-10Zr / HT-9 pin (single axial layer, full power from
t=0) twice:

  (a) COLD START — solver.set_hot_start(False) (the default). The pin
      begins at T0 = 293.15 K (room temperature); the real time integration
      itself resolves the thermal transient up to the hot steady state.

  (b) HOT START  — solver.set_hot_start(True). Before t=0 of the real run,
      the solver (own backward-Euler integrator, no IDA) repeatedly steps
      with a fixed t=0 (boundary conditions pinned) and geometrically
      growing dt, with afterAcceptedStep() (GRSIS/Zr redistribution/
      cladding wastage/burnup bookkeeping) never called during the march —
      i.e. irradiation physics off — until temperatures stop changing, then
      re-primes the real run's state from that converged point. t=0 in the
      hot-start run is already at the hot steady state.

TEND is kept short (300 s) on purpose: U-Pu-Zr metal fuel has very high
thermal conductivity, so thermal steady state is reached within ~30-40 s —
far too short for burnup/fission-gas release/Zr redistribution to become
significant, i.e. "hot steady state, before any serious irradiation
effects", per the task description. bup/fgrel are printed for both runs to
show they stay negligible over this window.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'build'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import _fred_m_na as fred   # compiled extension — full pybind11 API (set_hot_start, ...)

T0 = 293.15        # K — room temperature (cold start initial condition)
T_INLET = 643.0    # K — sodium coolant inlet
QQV = 3.67e9       # W/m3 — ~20 kW/m linear power, full power from t=0
TEND = 300.0       # s — reaches thermal steady state well before irradiation matters
DTOUT = 10.0       # s — fine enough to resolve the cold-start rise


def make_solver(hot_start: bool):
    g = fred.FuelRodGeometry()
    g.nf = 5        # fuel radial nodes
    g.nc = 3        # cladding radial nodes
    g.nz = 1        # single axial layer

    g.rfi0 = 0.0
    g.rfo0 = 2.35e-3
    g.rci0 = 2.45e-3    # 100 µm as-fab gap
    g.rco0 = 2.92e-3
    g.dz0  = [0.343]
    g.vgp  = 2.5e-6
    g.ruff = 5e-6
    g.rufc = 5e-6
    g.build()

    fuel = fred.UPuZr(pu_weight_frac=0.19, zr_weight_frac=0.10, reference_density=15700.0)
    clad = fred.HT9(reference_density=7750.0)

    s = fred.FredMNaSolver(g, fuel, clad)
    s.set_power_density_history([0.0, 1e6], [QQV, QQV])
    s.set_hot_start(hot_start)

    s.set_coolant_channel(
        dhyd=4.0e-3, xarea=6.0e-5, flowr=0.30,
        T_inlet_times=[0.0, 1e6], T_inlet_vals=[T_INLET, T_INLET],
    )
    s.set_coolant_pressure(0.1)   # MPa

    s.set_initial_temperature(T0)
    s.set_initial_gas_pressure(0.1)  # MPa

    s.set_enable_heat_conduction(True)
    s.set_enable_stress_strain(True)
    s.set_enable_zr_redistribution(True)
    s.set_enable_clad_wastage(True)
    return s


results = {}
for label, hot in [("cold", False), ("hot", True)]:
    print(f"\n=== FRED-M-Na {label.upper()} START ===")
    solver = make_solver(hot)
    solver.run(TEND, DTOUT, all_steps=False)
    times = np.asarray(solver.times())
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
ax.set_title('FRED-M-Na: cold start vs hot start\n(before serious irradiation effects)')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'hot_vs_cold_start.png'), dpi=120)
plt.close(fig)
print(f"\nSaved plots/hot_vs_cold_start.png")
