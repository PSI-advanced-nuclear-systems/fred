"""
FRED-ROD HDF5 output example — heat conduction + stress-strain.

Demonstrates:
  1. Running a FRED-ROD simulation with HDF5 output via output_file=
  2. Running again with all_steps=True to capture every IDA time step
  3. Reading results.h5 back with h5py / numpy and plotting

Physics:
  - Solid MOX pellet, T91 cladding, He gap
  - Ramp power 0 → 2e8 W/m3 over 50 s, then hold for 950 s
  - Heat conduction ON, stress-strain ON

Results units (FRED defaults, as stored in the HDF5 file):
  Temperature : K
  Stress      : MPa
  Radii       : m
  Time        : s
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     '..', '..', '..', 'build')
sys.path.insert(0, BUILD)

import fred_input as fred
import fred_output as fo

# ============================================================
# Geometry
# ============================================================
g = fred.FuelRodGeometry()
g.nf   = 5
g.nc   = 3
g.nz   = 1
g.rfi0 = 0.0         # solid pellet
g.rfo0 = 4.2e-3      # m
g.rci0 = 4.3e-3      # m  (100 µm as-fab gap)
g.rco0 = 4.9e-3      # m
g.dz0  = [0.10]
g.vgp  = 1.0e-6      # m3
g.ruff = 5.0e-6      # m
g.rufc = 5.0e-6      # m
g.build()

# ============================================================
# Materials
# ============================================================
fuel = fred.MOX(pu_content=0.15)         # low-Pu MOX (Philipponneau)
clad = fred.T91()                         # T91 ferritic-martensitic steel
gap  = fred.He()                          # pure He gap

# ============================================================
# Solver setup
# ============================================================
T0_K = 668.0   # K  initial / coolant temperature

def make_solver():
    s = fred.FredRodSolver(g, fuel, clad, gap)
    s.set_power_density_history([0.0, 50.0, 1000.0], [0.0, 2e8, 2e8])
    s.set_coolant_temperature([0.0, 1000.0], [T0_K, T0_K])
    s.set_coolant_pressure(0.1)
    s.set_initial_temperature(T0_K)
    s.set_initial_gas_pressure(0.1)
    s.set_enable_heat_conduction(True)
    s.set_enable_stress_strain(True)
    s.set_tolerances(1e-6, 1e-8)
    return s

# ============================================================
# Run 1 — dtout mode (every 50 s), write results.h5
# ============================================================
print("=== Run 1: dtout mode ===")
solver = make_solver()
solver.run(tend=1000.0, dtout=50.0, output_file='results.h5')
print(f"Written: results.h5  ({len(solver.time_points())} output steps)")

# ============================================================
# Run 2 — all_steps mode (every internal IDA step), write debug.h5
# ============================================================
print("\n=== Run 2: all_steps mode ===")
solver2 = make_solver()
solver2.run(tend=100.0, dtout=100.0, output_file='debug.h5', all_steps=True)
print(f"Written: debug.h5  ({len(solver2.time_points())} output steps — every IDA step)")

# ============================================================
# Read results.h5 and plot
# ============================================================
print("\n=== Reading results.h5 ===")
r = fo.read_results_h5('results.h5')

nz = r['nz']
nf = r['nf']
nc = r['nc']
times  = r['time']                  # [n_steps] s
T      = r['T']                     # [n_steps, nz, nf+nc] K
peak_T = r['peak_T_fuel']           # [n_steps] K
gap_w  = r['gap_width']             # [n_steps] m   (axial avg)
pfc    = r['pfc']                   # [n_steps] MPa
rfo    = r['rfo']                   # [n_steps] m
rci    = r['rci']                   # [n_steps] m
sigh   = r['sigh_outer']            # [n_steps] MPa
y_rst  = r['restart_y']             # [n_steps, neq]
yp_rst = r['restart_yp']            # [n_steps, neq]

print(f"  n_steps={r['n_steps']}, nz={nz}, nf={nf}, nc={nc}")
print(f"  restart_y shape: {y_rst.shape}")
print(f"\nFinal state (t={times[-1]:.0f} s):")
print(f"  T_fuel_centre = {T[-1, 0, 0] - 273.15:.1f} °C")
print(f"  T_fuel_outer  = {T[-1, 0, nf-1] - 273.15:.1f} °C")
print(f"  T_clad_inner  = {T[-1, 0, nf] - 273.15:.1f} °C")
print(f"  T_clad_outer  = {T[-1, 0, nf+nc-1] - 273.15:.1f} °C")
print(f"  gap width     = {gap_w[-1]*1e6:.1f} µm")
print(f"  pfc           = {pfc[-1]:.3f} MPa")
print(f"  sigh_outer    = {sigh[-1]:.3f} MPa")

# ============================================================
# Plots
# ============================================================
os.makedirs('plots', exist_ok=True)

# Temperature profiles
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ax = axes[0]
r_nodes = np.concatenate([
    np.linspace(g.rfi0, g.rfo0, nf),
    np.linspace(g.rci0, g.rco0, nc),
]) * 1e3  # mm
ax.plot(r_nodes[:nf], T[-1, 0, :nf] - 273.15, 'r-o', label='Fuel (final)')
ax.plot(r_nodes[nf:], T[-1, 0, nf:] - 273.15, 'b-o', label='Cladding (final)')
ax.plot(r_nodes[:nf], T[0, 0, :nf] - 273.15, 'r--', alpha=0.4, label='Fuel (t=0)')
ax.set_xlabel('Radius [mm]')
ax.set_ylabel('Temperature [°C]')
ax.set_title('Radial temperature profile')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(times, peak_T - 273.15, 'r-')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Peak fuel temperature [°C]')
ax.set_title('Fuel centreline temperature vs time')
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('plots/temperature.png', dpi=120)
plt.close(fig)
print("Saved plots/temperature.png")

# Mechanical results
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

axes[0].plot(times, gap_w * 1e6, 'k-')
axes[0].set_xlabel('Time [s]')
axes[0].set_ylabel('Gap width [µm]')
axes[0].set_title('Radial gap width (axial avg)')
axes[0].grid(True, alpha=0.3)

axes[1].plot(times, pfc, 'g-')
axes[1].set_xlabel('Time [s]')
axes[1].set_ylabel('Contact pressure [MPa]')
axes[1].set_title('Pellet-cladding contact pressure')
axes[1].grid(True, alpha=0.3)

axes[2].plot(times, sigh, 'b-')
axes[2].set_xlabel('Time [s]')
axes[2].set_ylabel('Hoop stress [MPa]')
axes[2].set_title('Outer cladding hoop stress')
axes[2].grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('plots/mechanical.png', dpi=120)
plt.close(fig)
print("Saved plots/mechanical.png")

# Comparison of dtout vs all_steps time resolution
r2 = fo.read_results_h5('debug.h5')
fig, ax = plt.subplots()
ax.plot(r2['time'], r2['peak_T_fuel'] - 273.15, 'b.', ms=3, label='all_steps (every IDA step)')
# overlay dtout
idx_lt100 = times <= 100.0
ax.plot(times[idx_lt100], peak_T[idx_lt100] - 273.15, 'r-', lw=2, label='dtout=50 s')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Peak fuel temperature [°C]')
ax.set_title('dtout vs all_steps time resolution (first 100 s)')
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig('plots/dtout_vs_allsteps.png', dpi=120)
plt.close(fig)
print("Saved plots/dtout_vs_allsteps.png")

print("\nDone.")
