"""
FRED-M-Na: OpenMP multithreading demo.

Runs the same simulation (20 axial layers) single-threaded and with
threads=2, 5, 10, then:
  1. Confirms the per-layer OpenMP parallelization does not change results
     (plots the temperature-profile difference vs. the single-threaded
     baseline — should be ~0 everywhere).
  2. Plots wall-clock run time and speedup vs. thread count.

Geometry uses more radial nodes than the basic FRED-M-Na examples (nf=20,
nc=10) and a longer irradiation time (40 days) than simple_irradiation's
10 days, so each layer's per-step Newton solve does enough work for
threading to matter.

Unlike FRED-ROD/FRED-OX (which still solve one global dense IDA Newton
system per step — see those apps' own multithread_demo for why that caps
their speedup), FRED-M-Na's one-step backward-Euler integrator solves each
axial layer's own local dense Newton system independently (layers only
couple through 3 rod-scalar globals, held fixed for the duration of one
outer Picard sweep — see FredMNaSolver::solveStepBackwardEuler's header
comment). Both the residual assembly AND each layer's own linear algebra
are parallelized here, so this demo is expected to show better scaling.
"""

import sys
import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'build')
sys.path.insert(0, BUILD)

import _fred_m_na as fred

# ============================================================
# Geometry — 20 axial layers, finer radial mesh than simple_irradiation
# ============================================================
NZ = 20
NF = 20
NC = 10

TEND = 40.0 * 86400.0    # 40 days [s]
DTOUT = 4.0 * 86400.0    # output every 4 days
THREAD_COUNTS = [1, 2, 5, 10]


def make_solver():
    g = fred.FuelRodGeometry()
    g.nf = NF
    g.nc = NC
    g.nz = NZ
    g.rfi0 = 0.0
    g.rfo0 = 2.35e-3
    g.rci0 = 2.45e-3
    g.rco0 = 2.92e-3
    g.dz0 = [0.343 / NZ] * NZ
    g.vgp = 2.5e-6
    g.ruff = 5e-6
    g.rufc = 5e-6
    g.build()

    fuel = fred.UPuZr(pu_weight_frac=0.19, zr_weight_frac=0.10, reference_density=15700.0)
    clad = fred.HT9(reference_density=7750.0)

    s = fred.FredMNaSolver(g, fuel, clad)
    qqv = 3.67e9
    s.set_power_density_history([0.0, TEND], [qqv, qqv])
    s.set_hot_start(True)
    s.set_coolant_channel(dhyd=4.0e-3, xarea=6.0e-5, flowr=0.30,
                           T_inlet_times=[0.0, TEND], T_inlet_vals=[643.0, 643.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(673.0)
    s.set_initial_gas_pressure(0.1)
    s.set_enable_heat_conduction(True)
    s.set_enable_stress_strain(True)
    s.set_enable_zr_redistribution(True)
    s.set_enable_clad_wastage(True)
    return s


# ============================================================
# Run at each thread count
# ============================================================
results = {}
for n in THREAD_COUNTS:
    print(f"=== Running FRED-M-Na with threads={n} ===")
    solver = make_solver()
    t0 = time.perf_counter()
    solver.run(TEND, DTOUT, all_steps=False, threads=n)
    wall = time.perf_counter() - t0
    T_flat = np.array(solver.temperatures())  # [n_steps, nz*(nf+nc)]
    T = T_flat.reshape(T_flat.shape[0], NZ, NF + NC)
    results[n] = {'wall': wall, 'T': T}
    print(f"    wall time = {wall:.2f} s")

baseline_T = results[1]['T']

print("\n=== Summary ===")
print(f"{'threads':>8} {'wall [s]':>10} {'speedup':>9} {'max |dT| vs threads=1| [K]':>28}")
for n in THREAD_COUNTS:
    wall = results[n]['wall']
    speedup = results[1]['wall'] / wall
    maxdiff = np.max(np.abs(results[n]['T'] - baseline_T))
    results[n]['maxdiff'] = maxdiff
    print(f"{n:>8} {wall:>10.2f} {speedup:>9.2f} {maxdiff:>28.3e}")

# ============================================================
# Plots
# ============================================================
os.makedirs('plots', exist_ok=True)

# ---- Plot 1: runtime and speedup vs thread count ----
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

walls = [results[n]['wall'] for n in THREAD_COUNTS]
speedups = [results[1]['wall'] / results[n]['wall'] for n in THREAD_COUNTS]

ax = axes[0]
ax.plot(THREAD_COUNTS, walls, 'o-', color='#2b6cb0', linewidth=2, markersize=7)
ax.set_xlabel('threads')
ax.set_ylabel('Wall-clock run time [s]')
ax.set_title(f'FRED-M-Na run time (nz={NZ}, nf={NF}, nc={NC})')
ax.grid(True, alpha=0.3)
ax.set_xticks(THREAD_COUNTS)

ax = axes[1]
ax.plot(THREAD_COUNTS, speedups, 'o-', color='#2f855a', linewidth=2, markersize=7, label='measured')
ax.plot(THREAD_COUNTS, THREAD_COUNTS, '--', color='#a0aec0', linewidth=1.5, label='ideal (linear)')
ax.set_xlabel('threads')
ax.set_ylabel('Speedup vs. threads=1')
ax.set_title('Speedup')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xticks(THREAD_COUNTS)

fig.tight_layout()
fig.savefig('plots/runtime_speedup.png', dpi=120)
plt.close(fig)
print("Saved plots/runtime_speedup.png")

# ---- Plot 2: temperature-profile difference vs threads=1, at final time ----
blues = plt.cm.Blues(np.linspace(0.4, 0.95, len(THREAD_COUNTS) - 1))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

r_nodes = np.concatenate([
    np.linspace(0.0, 2.35e-3, NF),
    np.linspace(2.45e-3, 2.92e-3, NC),
]) * 1e3  # mm

ax = axes[0]
ax.plot(r_nodes, baseline_T[-1, 0, :] - 273.15, 'k-', linewidth=2, label='threads=1 (reference)')
ax.set_xlabel('Radius [mm]')
ax.set_ylabel('Temperature [°C]')
ax.set_title('Final radial temperature profile')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes[1]
for color, n in zip(blues, THREAD_COUNTS[1:]):
    dT = results[n]['T'][-1, 0, :] - baseline_T[-1, 0, :]
    ax.plot(r_nodes, dT, 'o-', color=color, markersize=3, linewidth=1.2,
             label=f'threads={n}')
ax.axhline(0.0, color='#718096', linewidth=1, linestyle='--')
ax.set_xlabel('Radius [mm]')
ax.set_ylabel('ΔT vs. threads=1 [K]')
ax.set_title('Temperature difference (multithreaded − serial)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(-1e-6, 1e-6)

fig.tight_layout()
fig.savefig('plots/temperature_difference.png', dpi=120)
plt.close(fig)
print("Saved plots/temperature_difference.png")

print("\nDone.")
