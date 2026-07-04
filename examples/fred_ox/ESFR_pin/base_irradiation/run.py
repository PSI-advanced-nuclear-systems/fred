"""
FRED-OX ESFR-SMART pin — base irradiation 600 days.

Geometry: nz=20 × 0.05 m, nf=22, nc=5; MOX fuel, T91 cladding, He gap
Duration: 600 days (5.184e7 s); dtout = 86400 s (1 day per output)

Snapshot saved to snapshots/esfr_pin_base_frame1.snapshot for transient restart.

Legacy FRED reference at t=600d (ref_results/outfrd000000000600):
  gpres  = 0.59687 MPa
  fggen  = 428.871 cm3 STP
  fgrel  = 101.969 cm3 STP
  fgrel% = 23.776 %
  T_peak_fuel (node 13) = 1826.6 degC
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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SNAP_DIR   = os.path.join(SCRIPT_DIR, 'snapshots')
PLOT_DIR   = os.path.join(SCRIPT_DIR, 'plots')
os.makedirs(SNAP_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

# ============================================================
# Geometry (from fred.inp cards 000003, 100001, 100002, 100003)
# ============================================================
g = fred.FuelRodGeometry()
g.nz   = 20
g.nf   = 22       # radial fuel nodes
g.nc   = 5        # radial cladding nodes
g.rfi0 = 1.56e-3  # fuel inner radius [m]
g.rfo0 = 4.68e-3  # fuel outer radius [m]
g.rci0 = g.rfo0 + 1.55e-4   # rci0 = rfo0 + dgap = 4.835e-3 m
g.rco0 = 5.335e-3            # cladding outer radius [m]
g.dz0  = [0.05] * g.nz      # axial node height [m]; total = 1.0 m
g.vgp  = 7.19e-5             # plenum volume [m3]
g.ruff = 2.0e-6              # fuel surface roughness [m]
g.rufc = 2.0e-6              # clad surface roughness [m]
g.build()

print(f"Geometry: rfi={g.rfi0*1e3:.2f} mm  rfo={g.rfo0*1e3:.2f} mm  "
      f"rci={g.rci0*1e3:.3f} mm  rco={g.rco0*1e3:.3f} mm")
print(f"As-fab gap: {(g.rci0-g.rfo0)*1e6:.1f} µm   Active length: {g.nz*0.05:.2f} m")

# ============================================================
# Materials (from fred.inp cards 100001, 100002, 100003)
# ============================================================
mox  = fred.FredOxMOX(pu_content=0.1799, rof0=10542.0, sto0=1.97)
clad = fred.FredOxT91(reference_density=7790.0)
gap  = fred.FredOxGapMaterial()

# ============================================================
# Solver
# ============================================================
solver = fred.FredOxSolver(g, mox, clad, gap)

# ============================================================
# Boundary conditions (from fred.inp cards 200001, 200002, 200006)
# Time points: t=0, t=1h, t=1d, t=300d, t=600d
# ============================================================
times_bc = [0.0, 3600.0, 86400.0, 2.592e7, 5.184e7]  # s

# --- Power density per axial layer [W/m3] ---
# Table rows: one per time; 20 values per row (layer 0=bottom … layer 19=top)
# Transposed below to (n_layers × n_times) as required by the API.
_qqv_rows = [
    # t = 0 s
    [0.0]*20,
    # t = 3600 s (1 h)
    [0.0]*20,
    # t = 86400 s (1 day)
    [1.18638807e+07, 1.44915439e+07, 1.93168298e+07, 2.76509758e+07, 4.38064530e+07,
     4.21763051e+08, 4.87705005e+08, 5.48096493e+08, 6.00331209e+08, 6.41682184e+08,
     6.71606438e+08, 6.88475029e+08, 6.93410118e+08, 6.84360444e+08, 6.62918751e+08,
     6.29760747e+08, 5.85984460e+08, 5.33713546e+08, 4.79656829e+08, 4.42022254e+08],
    # t = 2.592e7 s (300 days)
    [3.57860290e+07, 3.95543130e+07, 4.93158942e+07, 6.58478378e+07, 8.85251130e+07,
     4.36978569e+08, 5.01593238e+08, 5.59233625e+08, 6.04976708e+08, 6.43347324e+08,
     6.69832702e+08, 6.84456974e+08, 6.86266909e+08, 6.76613924e+08, 6.54955038e+08,
     6.21507445e+08, 5.77393303e+08, 5.26473807e+08, 4.71415593e+08, 4.33949945e+08],
    # t = 5.184e7 s (600 days)
    [5.89471601e+07, 6.37808924e+07, 7.85620258e+07, 1.02407312e+08, 1.31063405e+08,
     4.48465621e+08, 5.10377454e+08, 5.65423602e+08, 6.06050603e+08, 6.40777217e+08,
     6.64632156e+08, 6.74900519e+08, 6.74779856e+08, 6.63051480e+08, 6.40885813e+08,
     6.06014404e+08, 5.61357282e+08, 5.12875164e+08, 4.58927044e+08, 4.21352799e+08],
]
# Transpose to [n_layers × n_times]
qqv_per_layer = [[row[j] for row in _qqv_rows] for j in range(g.nz)]
solver.set_power_density_history_per_layer(times_bc, qqv_per_layer)

# --- Cladding outer temperature per axial layer [K] ---
_Tco_rows = [
    # t = 0 s
    [668.0]*20,
    # t = 3600 s (1 h)
    [668.0]*20,
    # t = 86400 s (1 day)
    [668.50, 668.77, 669.16, 669.72, 670.66,
     682.50, 691.52, 701.51, 712.29, 723.66,
     735.42, 747.32, 759.16, 770.69, 781.70,
     791.99, 801.40, 809.81, 817.27, 824.31],
    # t = 2.592e7 s (300 days)
    [669.19, 669.88, 670.82, 672.11, 673.85,
     685.46, 694.57, 704.57, 715.20, 726.40,
     737.92, 749.54, 761.05, 772.26, 782.96,
     792.95, 802.06, 810.23, 817.41, 824.20],
    # t = 5.184e7 s (600 days)
    [669.85, 670.95, 672.42, 674.40, 676.90,
     688.27, 697.44, 707.46, 717.99, 729.03,
     740.36, 751.70, 762.92, 773.81, 784.20,
     793.85, 802.62, 810.54, 817.47, 824.01],
]
Tco_per_layer = [[row[j] for row in _Tco_rows] for j in range(g.nz)]
solver.set_coolant_temperature_per_layer(times_bc, Tco_per_layer)

# --- Plenum temperature [K] ---
solver.set_plenum_temperature_history(
    times_bc, [668.00, 668.00, 668.50, 669.19, 669.85])

# --- Coolant pressure [MPa] ---
solver.set_coolant_pressure(0.2)

# ============================================================
# Initial conditions and physics
# ============================================================
solver.set_initial_temperature(668.0)
solver.set_initial_gas_pressure(0.1)   # He fill at fabrication

solver.set_swelling_multiplier(0.8)    # FUEL_SWEL 0.8
solver.set_enable_heat_conduction(True)
solver.set_enable_stress_strain(True)
solver.set_tolerances(1.0e-5, 1.0e-7)

# ============================================================
# Snapshot: save final state for transient restart
# ============================================================
solver.set_snapshot_prefix(os.path.join(SNAP_DIR, 'esfr_pin_base'))
print(f"Snapshot will be saved to: {SNAP_DIR}/esfr_pin_base_frame1.snapshot")

# ============================================================
# Run (600 days, 1 day per output → 600+1 steps)
# ============================================================
TEND  = 5.184e7   # s (600 days)
DTOUT = 86400.0   # s (1 day)

print(f"\nRunning FRED-OX base irradiation: {TEND/86400:.0f} days, dtout={DTOUT/3600:.0f} h ...")
print("Estimated run time: 2–3 hours\n")

solver.run(tend=TEND, dtout=DTOUT)

# ============================================================
# Post-processing
# ============================================================
R_GAS = 8.314  # J/(mol·K)

times  = np.array(solver.time_points())   # s
gpres  = np.array(solver.gas_pressure())  # MPa
fggen  = np.array(solver.fg_generated())  # mol
fgrel  = np.array(solver.fg_released())   # mol
gap_w  = np.array(solver.gap_width())     # m (axial avg)
bup    = np.array(solver.burnup())        # MWd/kgU (axial avg)
T_peak = np.array(solver.peak_fuel_temperature())  # K

# Convert fission gas: mol → cm3 STP (0°C, 1 atm)
fggen_cm3 = fggen * R_GAS * 293.15 / 1.0e5 * 1.0e6
fgrel_cm3 = fgrel * R_GAS * 293.15 / 1.0e5 * 1.0e6
fgrel_pct = fgrel / np.where(fggen > 1e-30, fggen, 1e-30) * 100.0

print(f"\n=== Final state (t={times[-1]/86400:.0f} days) ===")
print(f"  gpres   = {gpres[-1]:.5f}  MPa")
print(f"  fggen   = {fggen_cm3[-1]:.3f}   cm3 STP")
print(f"  fgrel   = {fgrel_cm3[-1]:.3f}   cm3 STP")
print(f"  fgrel%  = {fgrel_pct[-1]:.3f}   %")
print(f"  T_peak  = {T_peak[-1]-273.15:.1f} degC")
print(f"  bup_avg = {bup[-1]:.3f}   MWd/kgU")
print(f"  gap_avg = {gap_w[-1]*1e6:.2f}  µm")

# Compare with legacy FRED reference
print("\n=== Comparison with legacy FRED (ref_results/outfrd000000000600) ===")
ref = {
    'gpres':      0.596873,
    'fggen_cm3':  428.871,
    'fgrel_cm3':  101.969,
    'fgrel_pct':  23.776,
    'T_peak_C':   1826.6,
}
print(f"  gpres:   FRED-OX={gpres[-1]:.5f}   Legacy={ref['gpres']:.5f}   MPa  "
      f"(Δ={gpres[-1]-ref['gpres']:+.5f})")
print(f"  fggen:   FRED-OX={fggen_cm3[-1]:.3f}   Legacy={ref['fggen_cm3']:.3f}   cm3  "
      f"(Δ={fggen_cm3[-1]-ref['fggen_cm3']:+.3f})")
print(f"  fgrel:   FRED-OX={fgrel_cm3[-1]:.3f}   Legacy={ref['fgrel_cm3']:.3f}   cm3  "
      f"(Δ={fgrel_cm3[-1]-ref['fgrel_cm3']:+.3f})")
print(f"  fgrel%:  FRED-OX={fgrel_pct[-1]:.3f}   Legacy={ref['fgrel_pct']:.3f}   %    "
      f"(Δ={fgrel_pct[-1]-ref['fgrel_pct']:+.3f})")
print(f"  T_peak:  FRED-OX={T_peak[-1]-273.15:.1f}   Legacy={ref['T_peak_C']:.1f}   degC  "
      f"(Δ={T_peak[-1]-273.15-ref['T_peak_C']:+.1f})")

# ============================================================
# Plots
# ============================================================
t_days = times / 86400.0

fig, ax = plt.subplots()
ax.plot(t_days, gpres, 'b-')
ax.axhline(ref['gpres'], color='r', linestyle='--', label='Legacy FRED (600d)')
ax.set_xlabel('Time [days]')
ax.set_ylabel('Gas pressure [MPa]')
ax.set_title('Internal gas pressure — 600-day base irradiation')
ax.legend()
fig.savefig(os.path.join(PLOT_DIR, 'gas_pressure.png'), dpi=150)
plt.close(fig)

fig, ax = plt.subplots()
ax.plot(t_days, T_peak - 273.15, 'r-')
ax.axhline(ref['T_peak_C'], color='k', linestyle='--', label='Legacy FRED (600d)')
ax.set_xlabel('Time [days]')
ax.set_ylabel('Peak fuel temperature [°C]')
ax.set_title('Peak fuel centreline temperature — 600-day base irradiation')
ax.legend()
fig.savefig(os.path.join(PLOT_DIR, 'peak_T_fuel.png'), dpi=150)
plt.close(fig)

fig, ax = plt.subplots()
ax.plot(t_days, fggen_cm3, 'b-', label='Generated')
ax.plot(t_days, fgrel_cm3, 'g-', label='Released')
ax.axhline(ref['fggen_cm3'], color='b', linestyle='--', alpha=0.5)
ax.axhline(ref['fgrel_cm3'], color='g', linestyle='--', alpha=0.5, label='Legacy FRED (600d)')
ax.set_xlabel('Time [days]')
ax.set_ylabel('Fission gas [cm3 STP]')
ax.set_title('Fission gas generation and release — 600-day base irradiation')
ax.legend()
fig.savefig(os.path.join(PLOT_DIR, 'fission_gas.png'), dpi=150)
plt.close(fig)

fig, ax = plt.subplots()
ax.plot(t_days, gap_w * 1e6, 'k-')
ax.set_xlabel('Time [days]')
ax.set_ylabel('Axial-avg gap width [µm]')
ax.set_title('Gap width — 600-day base irradiation')
fig.savefig(os.path.join(PLOT_DIR, 'gap_width.png'), dpi=150)
plt.close(fig)

fig, ax = plt.subplots()
ax.plot(t_days, bup, 'm-')
ax.set_xlabel('Time [days]')
ax.set_ylabel('Axial-avg burnup [MWd/kgU]')
ax.set_title('Burnup — 600-day base irradiation')
fig.savefig(os.path.join(PLOT_DIR, 'burnup.png'), dpi=150)
plt.close(fig)

print(f"\nPlots saved to {PLOT_DIR}/")
print("Snapshot saved for transient restart: snapshots/esfr_pin_base_frame1.snapshot")
print("\nDone.")
