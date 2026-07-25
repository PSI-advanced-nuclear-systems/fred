"""
FRED-M-Na: File to show fuel thermal conductivity burn-up profile based on the Timpano benchmark
Use only a single axial node 
U-Pu-Zr metallic fuel (13.25 wt% Pu, 10 wt% Zr) / HT-9 cladding
Sodium-cooled fast reactor, based on ESFR-SIMPLE geometry

Input parameters from ref_input.txt.
Compared against: legacy FRED-M and SAS4A at t=0, 45, 135, 365, 730, 1088, 2176 days.

Outputs:
  - fred_m_na_base_irradiation.h5 : full simulation results (HDF5, written by C++ at each step)
"""

import sys, os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BUILD_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'build'))
sys.path.insert(0, BUILD_DIR)
import _fred_m_na as fred


# ---------------------------------------------------------------------------
# Power: constant, axially uniform 1e+09 W/m3 (broadcast to all layers).
# For reference, the Timpano benchmark max power density is 1.0327E+09 W/m3.
# ---------------------------------------------------------------------------
power_times = [0.0, 2.0e8]
power_qqv_uniform = [1.0e9, 1.0e9]

# ---------------------------------------------------------------------------
# Geometry (from ref_input.txt)
# ---------------------------------------------------------------------------
NF   = 11
NC   = 5
NZ   = 24
DZ0  = 0.03125     # m
RFI  = 0.0         # m (solid fuel)
RFO  = 3.863e-3    # m
DGAP = 6.38e-4     # m
RCI  = RFO + DGAP
RCO  = 5.025e-3    # m
RUFF = 1.0e-6
RUFC = 1.0e-6
VGP  = 7.6375e-5   # m3

g = fred.FuelRodGeometry()
g.nf   = NF
g.nc   = NC
g.nz   = NZ
g.rfi0 = RFI
g.rfo0 = RFO
g.rci0 = RCI
g.rco0 = RCO
g.dz0  = [DZ0] * NZ
g.vgp  = VGP
g.ruff = RUFF
g.rufc = RUFC
g.build()

# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------
fuel = fred.UPuZr(
    pu_weight_frac    = 0.1325,
    zr_weight_frac    = 0.10,
    reference_density = 15800.0,
)
clad = fred.HT9(reference_density=7634.5)

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
solver = fred.FredMNaSolver(g, fuel, clad)
solver.set_grsis_data_mode(fred.GrsisDataMode.FEAST)
solver.set_sodium_mode(fred.SodiumMode.TDependent)
solver.set_conductivity_model(fred.ConductivityModel.EmpiricalBurnup)
solver.set_power_density_history(power_times, power_qqv_uniform)
solver.set_coolant_channel(
    dhyd          = 3.887e-3,
    xarea         = 3.301e-5,
    flowr         = 0.162,
    T_inlet_times = [0.0, 2.0e8],
    T_inlet_vals  = [633.0, 633.0],
    corr          = fred.HtcCorrelation.Subbotin,
)
solver.set_coolant_pressure(0.3856)
solver.set_initial_temperature(633.0)
solver.set_initial_gas_pressure(0.101)
# FRED-M-Na's one-step integrator (fixed-dt backward-Euler, "always accept",
# no SUNDIALS IDA — see FredMNaSolver.hpp) needs an explicit internal step, and taking
# that as one giant physics step would under-resolve the threshold-driven
# fanis "soft" contact transition (~day 100-130) and other fast-varying
# quantities (gpres, FGR) comparison. 
solver.set_step_size(86400.0)

# ---------------------------------------------------------------------------
# Run: full depletion 0 → 2176 days, with HDF5 output (C++ streaming)
# ---------------------------------------------------------------------------
# number of second in 1 day = 24*60*60 = 86400
D2S = 86400.0
TEND = 2176 * D2S
# DTOUT, 1 second to capture hot zero power, then 1 day for the rest of the run 
DTOUT = [1] + [D2S]*2176

THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
OUTPUT_H5  = os.path.join(THIS_DIR, "fred_m_na_base_irradiation.h5")

# Single 1-day step from cold start → thermal equilibrium with negligible burnup.
# Used for the "beginning of life" temperature comparison panel instead of the
# literal t=0 (cold, 633 K flat) or t=90 d (first main output, carries some burnup).

# Full irradiation run (resets to t=0 internally, independent of the above).
solver.set_hot_start(True)
solver.set_output_file(OUTPUT_H5)
solver.run(TEND, DTOUT)
print("Run complete.")

# ---------------------------------------------------------------------------
# Read results from HDF5 (C++ streamed the data step-by-step during the run)
# ---------------------------------------------------------------------------
print(f"Reading {os.path.basename(OUTPUT_H5)} ...")
with h5py.File(OUTPUT_H5, "r") as f:
    times_s   = f["time"][:]                  # [s], shape (n_steps,)
    T_flat    = f["thermal/T"][:]             # [K], shape (n_steps, nz*(nf+nc))
    k_fuel    = f["thermal/k_fuel"][:]        # [W/m-K], shape (n_steps, nz*nf)
    bup_atpct = f["burnup/bup_atpct"][:]      # [at%] rod-average burnup

times_d = times_s / D2S
stride  = NF + NC

# Fuel centerline (innermost node) temperature per axial layer, (n_steps, NZ)
T_fi_K   = T_flat[:, 0::stride]
# Peak fuel centerline temperature = max over axial layers at each step
T_peak_K = T_fi_K.max(axis=1)
# Axial layer where the end-of-life peak sits (used for the radial profiles)
j_peak   = int(T_fi_K[-1].argmax())

# ---------------------------------------------------------------------------
# Plot 1: peak fuel centerline temperature vs time, with a secondary x-axis
# in rod-average burnup [at%] (FRED-M-Na's internal burnup metric, streamed
# to HDF5 as burnup/bup_atpct).  Burnup grows monotonically at constant
# power, so days <-> at% is a well-defined interpolation.
# ---------------------------------------------------------------------------
d2b = lambda d: np.interp(d, times_d, bup_atpct)   # days  -> at%
b2d = lambda b: np.interp(b, bup_atpct, times_d)   # at%   -> days

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(times_d, T_peak_K - 273.15, '-', lw=1.5)
ax.set_xlabel('time [days]')
ax.set_ylabel('peak fuel centerline temperature [°C]')
ax.set_title('U-Pu-Zr, empirical piecewise burnup conductivity model (Karahan 2009)')
ax.grid(alpha=0.3)
secax = ax.secondary_xaxis('top', functions=(d2b, b2d))
secax.set_xlabel('rod-average burnup [at%]')
fig.tight_layout()
fig.savefig(os.path.join(THIS_DIR, 'peak_centerline_temperature.png'), dpi=140)
print("Plot written: peak_centerline_temperature.png")

# ---------------------------------------------------------------------------
# Plot 2: radial fuel temperature profile and fuel thermal conductivity at
# t = 1 s, mid-life (1088 d), and end of life (2176 d), at the peak axial
# layer.  k_fuel is drawn stepwise (steps-mid) to represent the per-node
# values the solver actually uses.
# ---------------------------------------------------------------------------
r_mm = np.linspace(RFI, RFO, NF) * 1e3   # as-fabricated node radii [mm]

fig, (ax_T, ax_k) = plt.subplots(1, 2, figsize=(12, 5))
for t_target, lbl in [(1.0,        't = 1 s (hot start, ~zero burnup)'),
                      (1088 * D2S, 't = 1088 d (mid-life)'),
                      (2176 * D2S, 't = 2176 d (end of life)')]:
    i = int(np.argmin(np.abs(times_s - t_target)))
    b = bup_atpct[i]
    T_prof = T_flat[i, j_peak*stride : j_peak*stride + NF] - 273.15
    k_prof = k_fuel[i, j_peak*NF : (j_peak+1)*NF]
    ax_T.plot(r_mm, T_prof, 'o-', label=f'{lbl}  [{b:.1f} at%]')
    ax_k.step(r_mm, k_prof, where='mid', label=f'{lbl}  [{b:.1f} at%]')
ax_T.set_xlabel('radius [mm]')
ax_T.set_ylabel('fuel temperature [°C]')
ax_T.set_title(f'Radial fuel temperature (axial layer {j_peak+1}/{NZ})')
ax_T.grid(alpha=0.3); ax_T.legend()
ax_k.set_xlabel('radius [mm]')
ax_k.set_ylabel('fuel thermal conductivity [W/m·K]')
ax_k.set_title(f'Fuel conductivity per node (axial layer {j_peak+1}/{NZ})')
ax_k.grid(alpha=0.3); ax_k.legend()
fig.tight_layout()
fig.savefig(os.path.join(THIS_DIR, 'radial_fuel_k_temp.png'), dpi=140)
print("Plot written: radial_fuel_k_temp.png")

# Console summary used by legacy_fred_m_na/compare.py
print("\nSummary (peak axial layer, centerline):")
for t_target, lbl in [(1.0, '1 s'), (1000 * D2S, '1000 d'), (2176 * D2S, '2176 d')]:
    i = int(np.argmin(np.abs(times_s - t_target)))
    k_prof = k_fuel[i, j_peak*NF : (j_peak+1)*NF]
    print(f"  t={lbl:8s}  T_cl={T_fi_K[i, j_peak]:8.2f} K  "
          f"k(center)={k_prof[0]:.3f} W/m-K  k(edge)={k_prof[-1]:.3f} W/m-K  "
          f"bup={bup_atpct[i]:.2f} at%")