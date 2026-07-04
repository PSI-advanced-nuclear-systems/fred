"""
Heat conduction + stress-strain — FRED-ROD benchmark.

Physics:  heat conduction ON,  stress-strain ON
Scenario: full power (1e8 W/m3) from t=0; coolant held at 700 K, pcool = 5 MPa.

Operator-split coupling: IDA advances temperatures; solveMechanicalIC updates
strains/stresses at each output time using the latest temperature distribution.
The gap in the thermal solve uses as-fabricated geometry (frozen coupling).

Hot start: full power is applied from t=0 directly (no manual power ramp).
Earlier versions of this example needed a slow 0->1e8 W/m3 ramp over the
first 50 s purely to keep IDACalcIC from failing on the initial jump from a
cold, uniform T0 to full power. set_hot_start(True) replaces that workaround:
before the real time integration starts, the solver marches the rod forward
in pseudo-time at its t=0 boundary conditions (irradiation physics off, not
that FRED-ROD has any) until temperatures/stresses stop changing, and uses
that converged state as t=0 for the real run — so t=0 in the printed table
and plots below is already at (or very near) steady state.

Expected steady-state behaviour:
  - Fuel centreline: ~837 K  (fuel ΔT = Q*R²/(4k) ≈ 44 K above outer fuel surface)
  - Cladding heats from 700 K to ~793 K (gap + cladding resistance adds ~93 K)
  - Thermal expansion of cladding (T0=700 K → 793 K) causes Δr_co to increase
  - External pressure (5 MPa) compresses clad; hoop stress becomes less compressive
    as temperature rises (thermal expansion partially offsets pressure effect)

Δr_co includes thermal expansion from T_REF (≈293 K) to the current cladding
temperature, so it is large and positive (~18–19 µm) even though the gap is open
and the pressure is compressive.

Results I/O: this example writes results.h5 via set_output_file() and reads
it straight back with fred_output.read_results_h5() — that HDF5 route is the
default, recommended way to get results out of FRED 2.0 (see fred_output.py's
module docstring for the file layout). The older direct-from-solver
accessors (solver.time_points()/.temperatures()/...) still work and are also
exercised once below, purely to plot them against the HDF5 read-back and
show the two agree to float64 round-off.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'build'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import _fred_rod as fred    # compiled extension — full pybind11 API (set_hot_start, set_output_file, ...)
import fred_output as fo    # shared HDF5 reader/writer (fred_platform/python/fred_output.py)

# ── Geometry ─────────────────────────────────────────────────────────────────
g = fred.FuelRodGeometry()
g.nf   = 3
g.nc   = 2
g.nz   = 1
g.rfi0 = 0.0        # solid pellet
g.rfo0 = 4.2e-3     # m
g.rci0 = 4.3e-3     # m   (100 µm as-fab gap)
g.rco0 = 4.9e-3     # m
g.dz0  = [0.10]     # m   (10 cm axial slice)
g.vgp  = 1.0e-6     # m3
g.ruff = 5.0e-6     # m
g.rufc = 5.0e-6     # m
g.build()

# ── Materials ─────────────────────────────────────────────────────────────────
fuel = fred.DummyFuelPellet()
clad = fred.DummyCladding()
gap  = fred.DummyGapMaterial()

# ── Solver ───────────────────────────────────────────────────────────────────
solver = fred.FredRodSolver(g, fuel, clad, gap)

solver.set_enable_heat_conduction(True)
solver.set_enable_stress_strain(True)

T0 = 700.0  # K  initial temperature
solver.set_initial_temperature(T0)
solver.set_coolant_temperature([0.0, 1e6], [T0, T0])
# Full power from t=0 — see module docstring: set_hot_start(True) below
# replaces the manual power ramp previously needed here for IC stability.
solver.set_power_density_history([0.0, 1e6], [1.0e8, 1.0e8])
solver.set_hot_start(True)

solver.set_coolant_pressure(5.0)        # MPa
solver.set_initial_gas_pressure(0.1)    # MPa  reference fill-gas pressure at T_REF=293.15 K
solver.set_tolerances(1e-6, 1e-8)

# ── Run ──────────────────────────────────────────────────────────────────────
tend  = 200.0   # s — match legacy: steady state reached well before t=200s
dtout = 10.0    # s
here    = os.path.dirname(os.path.abspath(__file__))
h5_path = os.path.join(here, 'results.h5')
solver.set_output_file(h5_path)
solver.run(tend, dtout)

# ── Results ──────────────────────────────────────────────────────────────────
# Method 1 — direct from the solver (in-memory, no file needed). This is the
# original way every FRED 2.0 example read results, and it still works —
# kept here as a comment for reference; it is exercised once below (as plain
# code, not literally commented out) only so it can be plotted against
# Method 2 to prove the two routes agree exactly.
#
#   times  = solver.time_points()
#   T_arr  = solver.temperatures()
#   sigh   = solver.clad_outer_hoop_stress()
#   rco    = solver.clad_outer_radius()
times_direct = solver.time_points()
T_direct     = solver.temperatures()
sigh_direct  = np.asarray(solver.clad_outer_hoop_stress())
rco_direct   = np.asarray(solver.clad_outer_radius())

# Method 2 (default / recommended) — read back from the HDF5 file written by
# set_output_file() above. This is the route newcomers should use: the same
# file also works with any non-Python HDF5 tool, and is what every other
# example in this repo now reads from by default.
r     = fo.read_results_h5(h5_path)
times = r['time']
# fred_output's native HDF5 stream stores /thermal/T flat as
# [n_steps, nz*(nf+nc)]; with nz=1 (this example) that is already
# [n_steps, nf+nc], matching solver.temperatures()'s shape exactly.
T_arr = r['T']
sigh  = r['sigh_outer']
rco   = r['rco']

nf, nc = g.nf, g.nc
rco0 = g.rco0

print(f"\n{'t [s]':>7}  {'T_ctr [K]':>10}  {'T_oci [K]':>10}  {'σ_θ [MPa]':>10}  {'Δr_co [µm]':>11}")
print("-" * 62)
for k, t in enumerate(times):
    T_row = T_arr[k]
    T_ctr = T_row[0]          # fuel centre node
    T_oci = T_row[nf]         # clad inner surface node
    dr_co = (rco[k] - rco0) * 1e6
    print(f"{t:7.2f}  {T_ctr:10.2f}  {T_oci:10.2f}  {sigh[k]:10.3f}  {dr_co:11.4f}")

print(f"\nSteady-state centreline temperature: {T_arr[-1][0]:.1f} K")
print(f"Outer clad radius change at steady state: {(rco[-1]-rco0)*1e6:.4f} µm")
print(f"  (includes thermal expansion from T_REF ≈ 293 K to clad temperature)")
print(f"Outer clad hoop stress at steady state:   {sigh[-1]:.3f} MPa")

# Old FRED legacy values at steady state (t=200s)
# Source: examples/primitive_cases/dummy_materials/heat_conduction_stress_strain/
# Slow power ramp 0→1e8 W/m3 over 50s; temperatures and stresses are at steady state.
# temperature (C) iz:1  569.381  558.356  525.281  438.354  426.850 → converted to K
Tfuel_legacy = [842.531, 831.506, 798.431]
Tclad_legacy = [711.504, 700.000]
sigh_outer_legacy = -28.4636   # MPa (sig h, outer clad node) — from legacy outfrd t=200s
rco_legacy        = 4.91881e-3 # m   (rco at t=200s) — from legacy outfrd t=200s

print("\nComparison with legacy FRED (steady-state values):")
print(f"{'Node':>6}  {'T [K]':>10}  {'Legacy FRED [K]':>16}  {'Material':>12}")
print("-" * 50)
for i in range(nf):
    T_node = T_arr[-1, i]
    print(f"{i:6d}  {T_node:10.3f}  {Tfuel_legacy[i]:16.3f}  {'Fuel':>12}")
for j in range(nc):
    T_node = T_arr[-1, nf + j]
    print(f"{nf+j:6d}  {T_node:10.3f}  {Tclad_legacy[j]:16.3f}  {'Cladding':>12}")

print(f"\n{'Quantity':>25}  {'New FRED':>12}  {'Legacy FRED':>12}")
print("-" * 52)
print(f"{'sigma_theta_outer [MPa]':>25}  {sigh[-1]:12.3f}  {sigh_outer_legacy:12.3f}")
print(f"{'Delta_rco [um]':>25}  {(rco[-1]-rco0)*1e6:12.4f}  {(rco_legacy-rco0)*1e6:12.4f}")

# ── Verify: Method 1 (direct) vs Method 2 (HDF5) are identical ─────────────
T_direct_ctr = T_direct[:, 0]
T_h5_ctr     = T_arr[:, 0]
max_abs_diff_T    = float(np.max(np.abs(T_direct_ctr - T_h5_ctr)))
max_abs_diff_sigh = float(np.max(np.abs(sigh_direct - sigh)))
max_abs_diff_rco  = float(np.max(np.abs(rco_direct - rco)))
print("\nMethod 1 (direct-from-solver) vs Method 2 (HDF5 read-back), should be ~0"
      " (only float64 write/read round-off):")
print(f"  max|ΔT_centreline| = {max_abs_diff_T:.3e} K")
print(f"  max|Δσ_θ|          = {max_abs_diff_sigh:.3e} MPa")
print(f"  max|Δr_co|         = {max_abs_diff_rco:.3e} m")
assert max_abs_diff_T < 1.0e-4 and max_abs_diff_sigh < 1.0e-4 and max_abs_diff_rco < 1.0e-9, \
    "direct-from-solver and HDF5 read-back results diverged by more than float64 round-off"

# ── Plots: overlay Method 1 and Method 2 to show they coincide exactly ─────
plots_dir = os.path.join(here, 'plots')
os.makedirs(plots_dir, exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

axes[0].plot(times_direct, T_direct[:, 0], 'b-', lw=4, alpha=0.35, label='Method 1: direct (solver)')
axes[0].plot(times, T_arr[:, 0], 'k--', label='Method 2: HDF5 read-back')
axes[0].set_xlabel('Time [s]'); axes[0].set_ylabel('Fuel centreline T [K]')
axes[0].set_title('Temperature'); axes[0].legend(fontsize=8); axes[0].grid(True, alpha=0.3)

axes[1].plot(times_direct, sigh_direct, 'b-', lw=4, alpha=0.35, label='Method 1: direct (solver)')
axes[1].plot(times, sigh, 'k--', label='Method 2: HDF5 read-back')
axes[1].set_xlabel('Time [s]'); axes[1].set_ylabel(r'$\sigma_\theta$ outer clad [MPa]')
axes[1].set_title('Hoop stress'); axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.3)

axes[2].plot(times_direct, (rco_direct - rco0) * 1e6, 'b-', lw=4, alpha=0.35, label='Method 1: direct (solver)')
axes[2].plot(times, (rco - rco0) * 1e6, 'k--', label='Method 2: HDF5 read-back')
axes[2].set_xlabel('Time [s]'); axes[2].set_ylabel(r'$\Delta r_{co}$ [µm]')
axes[2].set_title('Clad outer radius'); axes[2].legend(fontsize=8); axes[2].grid(True, alpha=0.3)

fig.suptitle('Direct-from-solver (Method 1) vs HDF5 read-back (Method 2) — identical by construction')
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'direct_vs_hdf5.png'), dpi=120)
plt.close(fig)
print(f"\nSaved plots/direct_vs_hdf5.png")
