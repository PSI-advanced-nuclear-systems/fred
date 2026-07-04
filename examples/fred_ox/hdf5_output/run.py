"""
FRED-OX HDF5 output example — 1-day base irradiation (ESFR-SMART pin).

Demonstrates:
  1. Running a FRED-OX simulation with HDF5 output via output_file=
  2. Running a short all_steps=True debug run to capture every IDA step
  3. Reading results.h5 back with h5py / numpy and plotting all key variables

Results units (FRED defaults, as stored in the HDF5 file):
  Temperature    : K
  Gas pressure   : MPa
  Fission gas    : mol  (convert to cm3 STP: x R_GAS x 293.15/1e5 x 1e6)
  Burnup         : MWd/kgU
  Gap width      : m
  Time           : s
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
# Geometry (from legacy fred.inp)
# ============================================================
g = fred.FuelRodGeometry()
g.nz  = 20
g.nf  = 22
g.nc  = 5
g.rfi0 = 1.56e-3
g.rfo0 = 4.68e-3
g.rci0 = g.rfo0 + 1.55e-4
g.rco0 = 5.335e-3
g.dz0  = [0.05] * g.nz
g.vgp  = 7.19e-5
g.ruff = 2.0e-6
g.rufc = 2.0e-6
g.build()

# ============================================================
# Materials
# ============================================================
PU_CONT = 0.1799
STO0    = 1.97
ROF0    = 10542.0

mox  = fred.FredOxMOX(pu_content=PU_CONT, rof0=ROF0, sto0=STO0)
clad = fred.FredOxT91(reference_density=7790.0)
gap  = fred.FredOxGapMaterial()

# ============================================================
# Per-layer power densities at t=86400 s (from fred.inp)
# ============================================================
qqv_at_tend = [
    1.18638807e+07, 1.44915439e+07, 1.93168298e+07, 2.76509758e+07,
    4.38064530e+07, 4.21763051e+08, 4.87705005e+08, 5.48096493e+08,
    6.00331209e+08, 6.41682184e+08, 6.71606438e+08, 6.88475029e+08,
    6.93410118e+08, 6.84360444e+08, 6.62918751e+08, 6.29760747e+08,
    5.85984460e+08, 5.33713546e+08, 4.79656829e+08, 4.42022254e+08,
]
T_co_at_tend = [
    668.50, 668.77, 669.16, 669.72, 670.66,
    682.50, 691.52, 701.51, 712.29, 723.66,
    735.42, 747.32, 759.16, 770.69, 781.70,
    791.99, 801.40, 809.81, 817.27, 824.31,
]

# ============================================================
# Solver setup (shared between runs)
# ============================================================
def make_solver():
    s = fred.FredOxSolver(g, mox, clad, gap)

    times_power = [0.0, 3600.0, 86400.0]
    qqv_per_layer = [[0.0, 0.0, qqv_at_tend[j]] for j in range(g.nz)]
    s.set_power_density_history_per_layer(times_power, qqv_per_layer)

    times_Tco = [0.0, 3600.0, 86400.0]
    Tco_per_layer = [[668.0, 668.0, T_co_at_tend[j]] for j in range(g.nz)]
    s.set_coolant_temperature_per_layer(times_Tco, Tco_per_layer)

    s.set_plenum_temperature_history([0.0, 3600.0, 86400.0], [668.0, 668.0, 668.5])
    s.set_coolant_pressure(0.2)
    s.set_initial_temperature(668.0)
    s.set_initial_gas_pressure(0.1)
    s.set_swelling_multiplier(0.8)
    s.set_enable_heat_conduction(True)
    s.set_enable_stress_strain(True)
    s.set_tolerances(1.0e-5, 1.0e-7)
    return s

# ============================================================
# Run 1 — dtout=1800 s (30 min), write results.h5
# ============================================================
print("=== Run 1: dtout=1800 s ===")
solver = make_solver()
solver.run(tend=86400.0, dtout=1800.0, output_file='results.h5')
print(f"Written: results.h5  ({len(solver.time_points())} output steps)")

# ============================================================
# Run 2 — short all_steps debug run (first 3600 s)
# ============================================================
print("\n=== Run 2: all_steps mode, first 7200 s ===")
solver2 = make_solver()
solver2.run(tend=7200.0, dtout=3600.0, output_file='debug.h5', all_steps=True)
print(f"Written: debug.h5  ({len(solver2.time_points())} output steps — every IDA step)")

# ============================================================
# Read results.h5
# ============================================================
print("\n=== Reading results.h5 ===")
r = fo.read_results_h5('results.h5')

nz = r['nz']
nf = r['nf']
nc = r['nc']
times  = r['time']
T      = r['T']            # [n_steps, nz, nf+nc] K
peak_T = r['peak_T_fuel']  # [n_steps] K
gap_w  = r['gap_width']    # [n_steps] m
gpres  = r['gpres']        # [n_steps] MPa
fggen  = r['fggen']        # [n_steps] mol
fgrel  = r['fgrel']        # [n_steps] mol
bup    = r['bup']          # [n_steps] MWd/kgU

R_GAS = 8.314  # J/(mol·K)
# Convert mol → cm3 STP (0°C, 1 atm)
fggen_cm3 = fggen * R_GAS * 293.15 / 1.0e5 * 1.0e6
fgrel_cm3 = fgrel * R_GAS * 293.15 / 1.0e5 * 1.0e6
fgrel_pct = fgrel / np.where(fggen > 1e-30, fggen, 1e-30) * 100.0

print(f"  n_steps={r['n_steps']}, nz={nz}, nf={nf}, nc={nc}")
print(f"  restart_y shape: {r['restart_y'].shape}")
print(f"\nFinal state (t={times[-1]:.0f} s):")
print(f"  gpres   = {gpres[-1]:.5f} MPa")
print(f"  fggen   = {fggen_cm3[-1]:.5f} cm3 STP")
print(f"  fgrel   = {fgrel_cm3[-1]:.5f} cm3 STP  ({fgrel_pct[-1]:.3f} %)")
print(f"  T_peak  = {peak_T[-1] - 273.15:.1f} °C")
print(f"  bup_avg = {bup[-1]:.5f} MWd/kgU")
print(f"  gap_avg = {gap_w[-1]*1e6:.1f} µm")

# Legacy FRED reference
print("\n=== Comparison with legacy FRED (t≈86400 s) ===")
print(f"  gpres: FRED-OX={gpres[-1]:.5f}  Legacy=0.24326  MPa")
print(f"  fggen: FRED-OX={fggen_cm3[-1]:.5f}  Legacy=0.33828  cm3 STP")
print(f"  fgrel: FRED-OX={fgrel_cm3[-1]:.5f}  Legacy=0.01344  cm3 STP")

# ============================================================
# Plots
# ============================================================
os.makedirs('plots', exist_ok=True)
t_h = times / 3600.0

# Gas pressure
fig, ax = plt.subplots()
ax.plot(t_h, gpres, 'b-')
ax.axhline(0.24326, color='r', linestyle='--', label='Legacy FRED (final)')
ax.set_xlabel('Time [h]')
ax.set_ylabel('Gas pressure [MPa]')
ax.set_title('Internal gas pressure — base irradiation 1d')
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig('plots/gas_pressure.png', dpi=120)
plt.close(fig)
print("Saved plots/gas_pressure.png")

# Peak fuel temperature
fig, ax = plt.subplots()
ax.plot(t_h, peak_T - 273.15, 'r-')
ax.set_xlabel('Time [h]')
ax.set_ylabel('Peak fuel temperature [°C]')
ax.set_title('Fuel centreline temperature — base irradiation 1d')
ax.grid(True, alpha=0.3)
fig.savefig('plots/peak_T_fuel.png', dpi=120)
plt.close(fig)
print("Saved plots/peak_T_fuel.png")

# Fission gas
fig, ax = plt.subplots()
ax.plot(t_h, fggen_cm3, 'b-', label='Generated')
ax.plot(t_h, fgrel_cm3, 'g-', label='Released')
ax.set_xlabel('Time [h]')
ax.set_ylabel('Fission gas [cm3 STP]')
ax.set_title('Fission gas — base irradiation 1d')
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig('plots/fission_gas.png', dpi=120)
plt.close(fig)
print("Saved plots/fission_gas.png")

# Burnup and gap width
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(t_h, bup, 'k-')
axes[0].set_xlabel('Time [h]')
axes[0].set_ylabel('Burnup [MWd/kgU]')
axes[0].set_title('Average burnup')
axes[0].grid(True, alpha=0.3)

axes[1].plot(t_h, gap_w * 1e6, 'k-')
axes[1].set_xlabel('Time [h]')
axes[1].set_ylabel('Gap width [µm]')
axes[1].set_title('Axial-average gap width')
axes[1].grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig('plots/burnup_gap.png', dpi=120)
plt.close(fig)
print("Saved plots/burnup_gap.png")

# Axial temperature profile at final time
fig, ax = plt.subplots()
z_centres = np.array([sum(g.dz0[:j+1]) - g.dz0[j]/2 for j in range(nz)])
T_fuel_centre = T[-1, :, 0]    # fuel centreline per layer
T_clad_outer  = T[-1, :, -1]   # outer clad surface per layer
ax.plot(z_centres * 100, T_fuel_centre - 273.15, 'r-o', label='Fuel centre')
ax.plot(z_centres * 100, T_clad_outer  - 273.15, 'b-o', label='Clad outer')
ax.set_xlabel('Axial position [cm]')
ax.set_ylabel('Temperature [°C]')
ax.set_title('Axial temperature profile at t=86400 s')
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig('plots/axial_temperature.png', dpi=120)
plt.close(fig)
print("Saved plots/axial_temperature.png")

# Read debug.h5 to show all_steps granularity
r2 = fo.read_results_h5('debug.h5')
fig, ax = plt.subplots()
ax.plot(r2['time'], r2['peak_T_fuel'] - 273.15, 'b.', ms=2,
        label=f"all_steps ({r2['n_steps']} points)")
ax.axvline(3600, color='r', linestyle='--', alpha=0.5, label='dtout=3600 s')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Peak fuel temperature [°C]')
ax.set_title('all_steps: IDA internal step resolution (first 7200 s)')
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig('plots/allsteps_detail.png', dpi=120)
plt.close(fig)
print("Saved plots/allsteps_detail.png")

print("\nDone.")
