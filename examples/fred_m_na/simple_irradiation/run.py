"""
FRED-M-Na: Simple irradiation example
U-Pu-Zr metallic fuel / HT-9 cladding / sodium-cooled fast reactor

This example simulates a 10-day base irradiation of a U-19Pu-10Zr fuel pin
at steady-state power, matching typical EBR-II conditions.

Physics included:
  - Heat conduction (thermal ODEs via IDA)
  - Burnup accumulation
  - Fission gas release (Karahan 2009 empirical model for U-Pu-Zr)
  - Fuel swelling (solid fission products, Ogata 1999)
  - Sodium infiltration (per node, post-processing)
  - Zirconium redistribution (diffusion model, post-processing)
  - Cladding wastage (lanthanide diffusion into HT-9, post-processing)
  - Stress-strain (thermo-elastic, Karahan 2009 moduli)

Geometry: typical EBR-II Mark-III pin
  Fuel pellet: rfi=0, rfo=2.35 mm, nf=5 nodes
  Cladding: rci=2.45 mm, rco=2.92 mm, nc=3 nodes
  Active height: 343 mm (one axial layer)
  Gas plenum volume: 2.5e-6 m3

Operating conditions:
  Linear power: ~20 kW/m (qqv ~ 3.67e9 W/m3)
  Coolant inlet temperature: 643 K (370 C, typical SFR inlet)
  Coolant pressure: 0.1 MPa (low pressure sodium)
  Initial gas pressure: 0.1 MPa (He backfill)
"""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Path setup ──────────────────────────────────────────────────────────────
BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import _fred_m_na as fred

# ── Geometry ─────────────────────────────────────────────────────────────────
g = fred.FuelRodGeometry()
g.nf = 5        # fuel radial nodes
g.nc = 3        # cladding radial nodes
g.nz = 1        # single axial layer for simplicity

# EBR-II Mark-III pin dimensions [m]
g.rfi0 = 0.0          # solid fuel (no central hole)
g.rfo0 = 2.35e-3      # outer fuel radius
g.rci0 = 2.45e-3      # inner cladding radius (100 µm as-fab gap)
g.rco0 = 2.92e-3      # outer cladding radius

g.dz0  = [0.343]      # active height [m]
g.vgp  = 2.5e-6       # plenum volume [m3]
g.ruff = 5e-6         # fuel roughness [m]
g.rufc = 5e-6         # cladding roughness [m]
g.build()

# ── Materials ────────────────────────────────────────────────────────────────
# U-19wt%Pu-10wt%Zr fuel (composition close to EBR-II driver)
fuel = fred.UPuZr(
    pu_weight_frac  = 0.19,   # 19 wt% Pu
    zr_weight_frac  = 0.10,   # 10 wt% Zr
    reference_density = 15700.0,  # ~75% TD for metallic fuel [kg/m3]
)

# HT-9 ferritic-martensitic steel cladding
clad = fred.HT9(reference_density=7750.0)

print(f"Fuel theoretical density   : {fuel.theoretical_density():.0f} kg/m3")
print(f"Fuel reference density     : {fuel.reference_density():.0f} kg/m3")
print(f"Fuel k at 800K             : {fuel.thermal_conductivity(800):.2f} W/mK")
print(f"Fuel cp at 800K            : {fuel.heat_capacity(800):.0f} J/kgK")
print(f"Fuel CTE at 800K           : {fuel.thermal_expansion_strain(800)*1e5:.2f} × 10⁻⁵")
print(f"HT-9 k at 700K             : {clad.thermal_conductivity(700):.2f} W/mK")
print(f"HT-9 E at 700K             : {clad.youngs_modulus(700)/1e3:.1f} GPa  [note: youngs_modulus returns MPa]")
print()

# ── Solver setup ─────────────────────────────────────────────────────────────
solver = fred.FredMNaSolver(g, fuel, clad)

# Power density [W/m3] — full power from t=0. set_hot_start(True) below
# marches the pin to thermal/mechanical steady state at this power (with
# burnup/fission-gas/Zr-redistribution/wastage physics off) before the real
# 10-day irradiation clock starts, so no manual power ramp is needed.
qqv_target = 3.67e9   # ~20 kW/m linear power for this geometry
tend = 10.0 * 86400   # 10 days [s]
dtout = 86400.0       # output every day

t_hist = [0.0, tend]
qqv_hist = [qqv_target, qqv_target]
solver.set_power_density_history(t_hist, qqv_hist)
solver.set_hot_start(True)

# Sodium coolant: subchannel BC with inlet at 643 K.
# Representative EBR-II subchannel geometry: dhyd=4mm, xarea=6e-5 m2, flowr=0.3 kg/s.
T_inlet = 643.0  # K inlet temperature (= T_cool from original example)
solver.set_coolant_channel(
    dhyd          = 4.0e-3,    # m  — hydraulic diameter
    xarea         = 6.0e-5,    # m2 — flow area
    flowr         = 0.30,      # kg/s
    T_inlet_times = [0.0, tend],
    T_inlet_vals  = [T_inlet, T_inlet],
)
solver.set_coolant_pressure(0.1)  # 0.1 MPa

# Initial conditions
T0 = 673.0  # start at ~400 C
solver.set_initial_temperature(T0)
solver.set_initial_gas_pressure(0.1)  # He backfill at 0.1 MPa

# Enable all physics
solver.set_enable_heat_conduction(True)
solver.set_enable_stress_strain(True)
solver.set_enable_zr_redistribution(True)
solver.set_enable_clad_wastage(True)


print(f"FRED-M-Na: running {tend/86400:.0f}-day irradiation...")
print(f"  nf={g.nf}, nc={g.nc}, nz={g.nz}, neq={solver.neq_total()}")
print(f"  Power density: {qqv_target/1e9:.2f} GW/m3  (linear power ~{qqv_target * 3.14159*(g.rfo0**2-g.rfi0**2)/1000:.1f} kW/m)")
print()

# ── Run ────────────────────────────────────────────────────────────────────
h5_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results.h5')
solver.set_output_file(h5_path)
solver.run(tend, dtout, all_steps=False)

# ── Extract results ──────────────────────────────────────────────────────────
# Default (recommended) route: read back from the HDF5 file just written —
# works for any post-processing tool, not just this Python session, and is
# crash-safe (each dtout step is flushed to disk as the run progresses).
#
# The direct-from-solver accessors below still work and are the simpler
# choice for a short interactive session; kept here as a comment for
# reference (see examples/fred_rod/heat_conduction_stress_strain/run.py for
# a side-by-side demonstration that the two routes agree):
#
#   times   = np.array(solver.times()) / 86400.0
#   gpres   = np.array(solver.gas_pressure())
#   fggen   = np.array(solver.fg_generated())
#   fgrel   = np.array(solver.fg_released())
#   gap     = np.array(solver.gap_width()) * 1e6
#   bup     = np.array(solver.burnup())
#   T_peak  = np.array(solver.peak_fuel_temperature())
#   xwast   = np.array(solver.clad_wastage()) * 1e6
import fred_output as fo
r       = fo.read_results_h5(h5_path)
times   = r['time'] / 86400.0        # days
gpres   = r['gpres']                 # MPa
fggen   = r['fggen']                 # mol
fgrel   = r['fgrel']                 # mol
gap     = r['gap_width'] * 1e6       # µm
bup     = r['bup']                   # MWd/kgU
T_peak  = r['peak_T_fuel']           # K
xwast   = r['xwast'] * 1e6 if 'xwast' in r else np.zeros_like(times)  # µm

print("=" * 65)
print(f"{'Quantity':<35} {'Initial':>12} {'Final':>12}")
print("=" * 65)
print(f"{'Gas pressure (MPa)':<35} {gpres[0]:>12.4f} {gpres[-1]:>12.4f}")
print(f"{'FGR generated (mol)':<35} {fggen[0]:>12.5f} {fggen[-1]:>12.5f}")
print(f"{'FGR released (mol)':<35} {fgrel[0]:>12.5f} {fgrel[-1]:>12.5f}")

fgrpp = 100.0 * fgrel[-1] / fggen[-1] if fggen[-1] > 0 else 0.0
print(f"{'FGR release fraction (%)':<35} {'—':>12} {fgrpp:>12.2f}")
print(f"{'Gap width (µm)':<35} {gap[0]:>12.2f} {gap[-1]:>12.2f}")
print(f"{'Average burnup (MWd/kgU)':<35} {bup[0]:>12.4f} {bup[-1]:>12.4f}")
print(f"{'Peak fuel temperature (K)':<35} {T_peak[0]:>12.1f} {T_peak[-1]:>12.1f}")
print(f"{'Peak fuel temperature (C)':<35} {T_peak[0]-273.15:>12.1f} {T_peak[-1]-273.15:>12.1f}")
print(f"{'Max clad wastage (µm)':<35} {xwast[0]:>12.4f} {xwast[-1]:>12.4f}")
print("=" * 65)

# ── Plots ──────────────────────────────────────────────────────────────────
os.makedirs(os.path.join(os.path.dirname(__file__), 'plots'), exist_ok=True)
plot_dir = os.path.join(os.path.dirname(__file__), 'plots')

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
fig.suptitle('FRED-M-Na: U-19Pu-10Zr / HT-9 — 10-day irradiation', fontsize=12)

ax = axes[0, 0]
ax.plot(times, T_peak - 273.15, 'r-o', markersize=4)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Peak fuel temperature (°C)')
ax.set_title('Peak Fuel Temperature')
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.plot(times, gpres, 'b-o', markersize=4)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Gas pressure (MPa)')
ax.set_title('Internal Gas Pressure')
ax.grid(True, alpha=0.3)

ax = axes[1, 0]
# Convert to cm3 STP (1 mol at STP = 22400 cm3)
fggen_cm3 = fggen * 22400.0
fgrel_cm3 = fgrel * 22400.0
ax.plot(times, fggen_cm3, 'g-', label='Generated')
ax.plot(times, fgrel_cm3, 'r-', label='Released')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Fission gas (cm³ STP)')
ax.set_title('Fission Gas')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
ax.plot(times, gap, 'k-o', markersize=4)
ax.set_xlabel('Time (days)')
ax.set_ylabel('Gap width (µm)')
ax.set_title('Fuel-Cladding Gap Width')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'irradiation_results.png'), dpi=120)
plt.close()
print(f"\nPlots saved to: {plot_dir}/irradiation_results.png")
