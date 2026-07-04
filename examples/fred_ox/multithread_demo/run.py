"""
FRED-OX: OpenMP multithreading demo.

Two independent, additive sources of threading are exercised together by
threads=N here:
  1. The per-axial-layer residual EVALUATION loop (hand-rolled OpenMP,
     with per-thread FredOxMOX/FredOxGapMaterial/HeatConduction clones).
  2. SUNDIALS' own in-built NVECTOR_OPENMP N_Vector, now used for IDA's
     y/yp/tolerance/id vectors whenever threads>1 — this threads the WRMS
     norm, dot-product, and linear-combination kernels IDA itself runs on
     every step, independent of what the per-layer residual loop computes.

Runs the same simulation (20 axial layers) single-threaded and with
threads=2, 5, 10, then:
  1. Confirms neither source of parallelism changes results beyond
     floating-point reduction-order noise (plots the temperature-profile
     difference vs. the single-threaded baseline — should be ~0, i.e.
     machine-precision-scale, everywhere).
  2. Plots wall-clock run time and speedup vs. thread count.

Geometry uses more radial nodes than the basic FRED-OX examples (nf=20,
nc=8) so each layer's residual evaluation does enough work for threading
to matter.

Uses the raw compiled module (`_fred_ox`) directly rather than
python/fred_ox.py's OO wrapper: that wrapper's `FredOxMOX` class has a
pre-existing, threading-unrelated unit bug (`reference_density_frac` is
documented as a fraction of theoretical density but passed straight
through to the C++ constructor's `rof0` argument, which expects an
absolute density in kg/m3) — bypassed here by supplying the absolute
density directly, matching how the raw module's constructor is documented.
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

import _fred_ox as fred_ox

# ============================================================
# Geometry — 20 axial layers, finer radial mesh than the basic examples
# ============================================================
NZ = 20
NF = 20
NC = 8


def make_solver():
    g = fred_ox.FuelRodGeometry()
    g.nf = NF
    g.nc = NC
    g.nz = NZ
    g.rfi0 = 0.0
    g.rfo0 = 4.2e-3
    g.rci0 = 4.3e-3
    g.rco0 = 4.9e-3
    g.dz0 = [0.10] * NZ
    g.vgp = 1.0e-6
    g.ruff = 5.0e-6
    g.rufc = 5.0e-6
    g.build()

    fuel = fred_ox.FredOxMOX(0.18, 10400.0, 1.97)  # pu_content, rof0 [kg/m3], sto0
    clad = fred_ox.FredOxAIM1()
    gap = fred_ox.FredOxGapMaterial()

    s = fred_ox.FredOxSolver(g, fuel, clad, gap)
    s.set_power_density_history([0.0, 50.0, 2000.0], [0.0, 2e8, 2e8])
    s.set_coolant_temperature([0.0, 2000.0], [568.0, 568.0])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(568.0)
    s.set_initial_gas_pressure(0.1)
    s.set_tolerances(1e-7, 1e-9)
    return s


TEND = 2000.0
DTOUT = 200.0
THREAD_COUNTS = [1, 2, 5, 10]

# ============================================================
# Run at each thread count
# ============================================================
results = {}
for n in THREAD_COUNTS:
    print(f"=== Running FRED-OX with threads={n} ===")
    solver = make_solver()
    t0 = time.perf_counter()
    solver.run(TEND, DTOUT, False, n)
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

print(
    "\nNote: like FRED-ROD, FRED-OX still solves ONE global dense IDA "
    "Newton system per step via SUNDIALS' dense linear solver "
    "(SUNLinSol_Dense, O(neq^3)) — SUNDIALS has no in-built threaded dense "
    "direct solver, so that factorization stays serial regardless of "
    "threads=. What threads=N does parallelize: (a) the O(nz) independent "
    "per-layer residual evaluations (hand-rolled OpenMP) and (b) IDA's own "
    "O(neq) vector kernels (WRMS norms, dot products, linear combinations) "
    "via SUNDIALS' in-built NVECTOR_OPENMP N_Vector. Both are modest "
    "fractions of total wall time next to the dense factorization/solve "
    "they can't touch, so speedup is modest relative to FRED-M-Na (see "
    "examples/fred_m_na/multithread_demo), whose per-layer Newton solve "
    "keeps each layer's own linear algebra local, not funneled through a "
    "single shared system."
)

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
ax.set_title(f'FRED-OX run time (nz={NZ}, nf={NF}, nc={NC})')
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
    np.linspace(0.0, 4.2e-3, NF),
    np.linspace(4.3e-3, 4.9e-3, NC),
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
# Difference is machine-precision-scale (~1e-12 K, not exactly 0): NVECTOR_
# OPENMP's WRMS-norm/dot-product kernels sum in a thread-count-dependent
# order, a different but equally valid floating-point rounding, not a
# correctness issue. Fix the y-range so "flat line at ~0" is visibly a
# real, resolved near-zero rather than an accident of autoscaling.
ax.set_ylim(-1e-6, 1e-6)

fig.tight_layout()
fig.savefig('plots/temperature_difference.png', dpi=120)
plt.close(fig)
print("Saved plots/temperature_difference.png")

print("\nDone.")
