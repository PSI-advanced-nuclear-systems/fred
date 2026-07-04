"""
FRED-OX base irradiation example — 1-day irradiation (ESFR-SMART pin)

Key parameters:
  fmat = mox     pucont = 0.1799   fden = 10542 kg/m3
  gmat = he      pin = 0.1 MPa     vpl = 7.19e-5 m3
  cmat = t91     rco = 5.335e-3 m  roc0 = 7790 kg/m3
  nz = 20  dz = 0.05 m  (total active length = 1.0 m)
  tend = 86400 s (1 day)  dtout = 1800 s
  Options: FGR, FUEL_CREEP (C1=1e-36), FUEL_SWEL 0.8, CLAD_CREEP, CLAD_PLAS
  Note: CLAD_CREEP and CLAD_PLAS have zero rate for T91 in legacy FRED.
        FUEL_CREEP has near-zero rate for C1=1e-36.
        Dominant effects: swelling (0.8*MATPRO), FGR, burnup-dependent conductivity.
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'build'))
import fred_ox as fred

# ============================================================
# Geometry (from legacy fred.inp)
# ============================================================
g = fred.FuelRodGeometry()
g.nz  = 20
g.nf  = 22      # fuel radial nodes (legacy: nf=22)
g.nc  = 5       # cladding radial nodes (legacy: nc=5)
g.rfi0 = 1.56e-3   # inner fuel radius [m]
g.rfo0 = 4.68e-3   # outer fuel radius [m]
g.rci0 = g.rfo0 + 1.55e-4  # rci0 = rfo0 + dgap = 4.68e-3 + 1.55e-4 m
g.rco0 = 5.335e-3  # outer clad radius [m]
g.dz0  = [0.05] * g.nz  # axial layer height [m]
g.vgp  = 7.19e-5   # gas plenum volume [m3]
g.ruff = 2.0e-6    # fuel surface roughness [m]
g.rufc = 2.0e-6    # clad surface roughness [m]
g.build()

print(f"Geometry: rfi={g.rfi0*1e3:.2f} rfo={g.rfo0*1e3:.2f} rci={g.rci0*1e3:.3f} rco={g.rco0*1e3:.3f} [mm]")
print(f"Gap as-fabricated: {(g.rci0-g.rfo0)*1e6:.1f} µm")

# ============================================================
# Materials
# ============================================================
PU_CONT = 0.1799
STO0    = 1.97
ROF0    = 10542.0  # kg/m3

mox = fred.FredOxMOX(pu_content=PU_CONT, rof0=ROF0, sto0=STO0)
clad = fred.FredOxT91(reference_density=7790.0)
gap  = fred.FredOxGapMaterial()

print(f"MOX melting temp: {mox.melting_temperature():.1f} K")
print(f"MOX k(700K): {mox.thermal_conductivity(700):.3f} W/mK")

# ============================================================
# Solver setup
# ============================================================
solver = fred.FredOxSolver(g, mox, clad, gap)

# Power off during construction, then ramp to full power over 1 hour
# (legacy FRED: power=0 for t=0..3600s, then ramp to full at t=86400s)
# Per-layer power densities at t=86400s (from fred.inp card 200002):
qqv_at_tend = [
    1.18638807e+07,  # layer  0 (z=0.025 m)
    1.44915439e+07,  # layer  1
    1.93168298e+07,  # layer  2
    2.76509758e+07,  # layer  3
    4.38064530e+07,  # layer  4
    4.21763051e+08,  # layer  5
    4.87705005e+08,  # layer  6
    5.48096493e+08,  # layer  7
    6.00331209e+08,  # layer  8
    6.41682184e+08,  # layer  9
    6.71606438e+08,  # layer 10
    6.88475029e+08,  # layer 11
    6.93410118e+08,  # layer 12 (peak)
    6.84360444e+08,  # layer 13
    6.62918751e+08,  # layer 14
    6.29760747e+08,  # layer 15
    5.85984460e+08,  # layer 16
    5.33713546e+08,  # layer 17
    4.79656829e+08,  # layer 18
    4.42022254e+08,  # layer 19
]

times_power = [0.0, 3600.0, 86400.0]  # s
# Per-layer power history: zero until t=3600s, then ramp to full
qqv_per_layer = []
for j in range(g.nz):
    qqv_per_layer.append([0.0, 0.0, qqv_at_tend[j]])

solver.set_power_density_history_per_layer(times_power, qqv_per_layer)

# Outer cladding temperature (per-layer, from fred.inp card 200001)
T_co_at_t0 = [668.0] * 20  # uniform at t=0 and t=1h
T_co_at_tend = [
    668.50, 668.77, 669.16, 669.72, 670.66,
    682.50, 691.52, 701.51, 712.29, 723.66,
    735.42, 747.32, 759.16, 770.69, 781.70,
    791.99, 801.40, 809.81, 817.27, 824.31,
]
times_Tco = [0.0, 3600.0, 86400.0]
Tco_per_layer = []
for j in range(g.nz):
    Tco_per_layer.append([T_co_at_t0[j], T_co_at_t0[j], T_co_at_tend[j]])

solver.set_coolant_temperature_per_layer(times_Tco, Tco_per_layer)

# Gas plenum temperature (from fred.inp card 200006)
solver.set_plenum_temperature_history([0.0, 3600.0, 86400.0], [668.0, 668.0, 668.5])

# Coolant pressure
solver.set_coolant_pressure(0.2)  # MPa (from card 200005)

# Initial conditions
T0_K = 668.0
solver.set_initial_temperature(T0_K)
solver.set_initial_gas_pressure(0.1)  # MPa (from card 100002: pin=0.1)

# Physics options
solver.set_swelling_multiplier(0.8)  # FUEL_SWEL 0.8

solver.set_enable_heat_conduction(True)
solver.set_enable_stress_strain(True)

solver.set_tolerances(1.0e-5, 1.0e-7)

# ============================================================
# Run
# ============================================================
TEND  = 86400.0  # s  (1 day)
DTOUT = 1800.0   # s  (30 min)

print(f"\nRunning FRED-OX: tend={TEND:.0f} s, dtout={DTOUT:.0f} s ...")
solver.run(tend=TEND, dtout=DTOUT)

# ============================================================
# Post-processing
# ============================================================
times = np.array(solver.time_points())
T_all = solver.temperatures()  # shape [n_times, nz*(nf+nc)]
gpres = np.array(solver.gas_pressure())
fggen = np.array(solver.fg_generated())
fgrel = np.array(solver.fg_released())
gap_w = np.array(solver.gap_width())
bup   = np.array(solver.burnup())
T_peak= np.array(solver.peak_fuel_temperature())

nf, nc, nz = g.nf, g.nc, g.nz

# Convert fission gas to cm3 STP (mol -> cm3 at STP: 1 mol = 22400 cm3 at 0°C, 1 atm)
R_GAS = 8.314  # J/(mol·K)
fggen_cm3 = fggen * R_GAS * 293.15 / 1.0e5 * 1.0e6  # mol -> cm3 STP
fgrel_cm3 = fgrel * R_GAS * 293.15 / 1.0e5 * 1.0e6
fgrel_pct = fgrel / np.where(fggen > 1e-30, fggen, 1e-30) * 100.0

print(f"\n=== Final state (t={times[-1]:.0f} s) ===")
print(f"  gpres   = {gpres[-1]:.5e} MPa")
print(f"  fggen   = {fggen_cm3[-1]:.5e} cm3")
print(f"  fgrel   = {fgrel_cm3[-1]:.5e} cm3")
print(f"  fgrel%  = {fgrel_pct[-1]:.5e} %")
print(f"  T_peak  = {T_peak[-1]-273.15:.1f} °C")
print(f"  gap_avg = {gap_w[-1]*1e6:.1f} µm")
print(f"  bup_avg = {bup[-1]:.5e} MWd/kgU")

# Compare with legacy FRED reference at t=86400s:
print("\n=== Comparison with legacy FRED (irradiate_10d/outfrd000000000048) ===")
legacy_gpres = 0.24326    # MPa
legacy_fggen = 0.33828    # cm3 STP
legacy_fgrel = 0.013437   # cm3 STP
legacy_fgrel_pct = 3.972  # %
print(f"  gpres:  FRED-OX={gpres[-1]:.5f}  Legacy={legacy_gpres:.5f}  MPa")
print(f"  fggen:  FRED-OX={fggen_cm3[-1]:.5f}  Legacy={legacy_fggen:.5f}  cm3")
print(f"  fgrel:  FRED-OX={fgrel_cm3[-1]:.5f}  Legacy={legacy_fgrel:.5f}  cm3")
print(f"  fgrel%: FRED-OX={fgrel_pct[-1]:.3f}  Legacy={legacy_fgrel_pct:.3f}  %")

# ============================================================
# Plots
# ============================================================
os.makedirs('plots', exist_ok=True)

# Gas pressure vs time
fig, ax = plt.subplots()
ax.plot(times/3600, gpres, 'b-', label='FRED-OX')
ax.axhline(legacy_gpres, color='r', linestyle='--', label='Legacy FRED (final)')
ax.set_xlabel('Time [h]')
ax.set_ylabel('Gas pressure [MPa]')
ax.set_title('Internal gas pressure — base irradiation 1d')
ax.legend()
fig.savefig('plots/gas_pressure.png', dpi=120)
plt.close(fig)
print("Saved plots/gas_pressure.png")

# Peak fuel temperature vs time
fig, ax = plt.subplots()
ax.plot(times/3600, T_peak - 273.15, 'r-')
ax.set_xlabel('Time [h]')
ax.set_ylabel('Peak fuel temperature [°C]')
ax.set_title('Peak fuel centreline temperature — base irradiation 1d')
fig.savefig('plots/peak_T_fuel.png', dpi=120)
plt.close(fig)
print("Saved plots/peak_T_fuel.png")

# Fission gas generation and release vs time
fig, ax = plt.subplots()
ax.plot(times/3600, fggen_cm3, 'b-', label='Generated')
ax.plot(times/3600, fgrel_cm3, 'g-', label='Released')
ax.set_xlabel('Time [h]')
ax.set_ylabel('Fission gas [cm3 STP]')
ax.set_title('Fission gas — base irradiation 1d')
ax.legend()
fig.savefig('plots/fission_gas.png', dpi=120)
plt.close(fig)
print("Saved plots/fission_gas.png")

# Gap width vs time
fig, ax = plt.subplots()
ax.plot(times/3600, gap_w * 1e6, 'k-')
ax.set_xlabel('Time [h]')
ax.set_ylabel('Axial-average gap width [µm]')
ax.set_title('Gap width — base irradiation 1d')
fig.savefig('plots/gap_width.png', dpi=120)
plt.close(fig)
print("Saved plots/gap_width.png")

print("\nDone.")
