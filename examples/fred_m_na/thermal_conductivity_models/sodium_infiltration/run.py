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
solver.set_conductivity_model(fred.ConductivityModel.DetailedNaSodium)
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
    sw_close  = f["grsis/swclose"][:]         # [-] closed-bubble swelling per node
    sw_open   = f["grsis/swopen"][:]          # [-] open-bubble swelling per node
    zr_wf_n   = f["burnup/zr_wf"][:]          # [-] local Zr weight fraction per node

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
ax.set_title('U-Pu-Zr, detailed Na-infiltration conductivity model')
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

# ---------------------------------------------------------------------------
# Plot 3: radial irradiation scaling factor  k / k_fresh(T_local).
#
# Dividing the solver's per-node k by the fresh-fuel correlation evaluated at
# each node's own temperature (nominal composition) strips out the k(T) trend
# and isolates the irradiation knock-down the detailed model applies:
#
#   - inner core (r < 0.6 r_fo, sodium gate closed, f_inf = 0):
#     factor = (1 - p_eff)^1.5 x (Zr-redistribution composition effect).
#     Fission-gas swelling (per-node GRSIS) and gamma-phase Zr enrichment
#     both peak at the hot centerline, so the factor RISES outward --
#     which is why k increases from the center to ~2.5 mm despite the
#     falling temperature.
#   - outer band (r > 0.6 r_fo, alpha/alpha+delta): the sodium gate opens;
#     Na-logged pores conduct (Maxwell-Eucken P_c) and partially recover
#     the factor, giving the step structure.
#
# Note: k_fresh here uses the NOMINAL composition, so the plotted factor
# bundles the porosity term with the local Zr-redistribution effect on
# k_fresh (per-node Zr is not streamed to HDF5).
# ---------------------------------------------------------------------------
def kfresh_nominal(T_K, x_pu=0.1325, x_zr=0.10):
    A = 17.5 * ((1.0 - 2.23 * x_zr) / (1.0 + 1.61 * x_zr) - 2.62 * x_pu)
    B = 1.54e-2 * (1.0 + 0.061 * x_zr) / (1.0 + 1.61 * x_zr)
    C = 9.38e-6 * (1.0 - 2.7 * x_pu)
    return A + B * T_K + C * T_K**2

fig, (ax, ax_d) = plt.subplots(1, 2, figsize=(13, 5))
for t_target, lbl in [(1.0,        't = 1 s (hot start, ~zero burnup)'),
                      (1088 * D2S, 't = 1088 d (mid-life)'),
                      (2176 * D2S, 't = 2176 d (end of life)')]:
    i = int(np.argmin(np.abs(times_s - t_target)))
    b = bup_atpct[i]
    T_prof = T_flat[i, j_peak*stride : j_peak*stride + NF]
    k_prof = k_fuel[i, j_peak*NF : (j_peak+1)*NF]
    factor = k_prof / kfresh_nominal(T_prof)
    ax.step(r_mm, factor, where='mid', label=f'{lbl}  [{b:.1f} at%]')
ax.axvline(0.6 * RFO * 1e3, color='k', ls=':', lw=1.0, alpha=0.7)
ax.text(0.6 * RFO * 1e3 + 0.05, 0.15, 'radius gate $0.6\\,r_{fo}$',
        rotation=90, fontsize=9, va='bottom')
# Sodium infiltration needs BOTH gates open: radius (r > 0.6 r_fo) AND phase
# (alpha or alpha+delta, i.e. T below the alpha/beta boundary "line3" of the
# Fphase diagram).  At this composition line3 ~ 888 K, which the EOL profile
# only crosses well outside 0.6 r_fo -- between the two markers the fuel is
# still beta/gamma, f_inf = 0, and there is no sodium recovery.
pu100 = 0.1325 * 100.0 / 19.0
T_line3 = (1.0 - pu100) * 935.15 + pu100 * 868.15
i_eol = int(np.argmin(np.abs(times_s - 2176 * D2S)))
T_eol = T_flat[i_eol, j_peak*stride : j_peak*stride + NF]
r_phase_mm = np.interp(T_line3, T_eol[::-1], r_mm[::-1])
ax.axvline(r_phase_mm, color='g', ls='--', lw=1.0, alpha=0.7)
ax.text(r_phase_mm + 0.05, 0.15,
        f'phase gate $T<{T_line3:.0f}$ K (EOL)',
        rotation=90, fontsize=9, va='bottom', color='g')
ax.set_xlabel('radius [mm]')
ax.set_ylabel('irradiation scaling factor  $k / k_{fresh}(T)$  [-]')
ax.set_ylim(0.0, 1.1)
ax.set_title(f'Irradiation knock-down (axial layer {j_peak+1}/{NZ})')
ax.grid(alpha=0.3); ax.legend(loc='lower left')

# Right panel: multiplicative decomposition of the EOL factor.
#   factor = porosity term x Zr-composition term x sodium term
#   - porosity   : (1 - sw_gas)^1.5 with sw_gas = swclose + swopen (GRSIS)
#   - composition: k_fresh(T, local Zr) / k_fresh(T, nominal Zr)
#   - sodium     : the residual factor/(por*comp); = 1 where f_inf = 0,
#                  > 1 where Na-logged pores conduct (P_c and reduced p_eff)
sc_eol  = sw_close[i_eol].reshape(NZ, NF)[j_peak]
so_eol  = sw_open [i_eol].reshape(NZ, NF)[j_peak]
zr_eol  = zr_wf_n [i_eol].reshape(NZ, NF)[j_peak]
k_eol   = k_fuel  [i_eol, j_peak*NF:(j_peak+1)*NF]
fac_eol = k_eol / kfresh_nominal(T_eol)
por_t   = (1.0 - (sc_eol + so_eol))**1.5
comp_t  = kfresh_nominal(T_eol, x_zr=zr_eol) / kfresh_nominal(T_eol)
na_t    = fac_eol / (por_t * comp_t)

ax_d.step(r_mm, fac_eol, where='mid', color='k', lw=2, label='total factor')
ax_d.step(r_mm, por_t,   where='mid', label='porosity  $(1-sw_{gas})^{1.5}$')
ax_d.step(r_mm, comp_t,  where='mid', label='Zr redistribution  $k_f(Zr_{loc})/k_f(Zr_{nom})$')
ax_d.step(r_mm, na_t,    where='mid', label='sodium recovery  ($P_c$, reduced $p_{eff}$)')
ax_d.axhline(1.0, color='gray', lw=0.6)
ax_d.axvline(0.6 * RFO * 1e3, color='k', ls=':', lw=1.0, alpha=0.7)
ax_d.axvline(r_phase_mm, color='g', ls='--', lw=1.0, alpha=0.7)
ax_d.set_xlabel('radius [mm]')
ax_d.set_ylabel('multiplicative contribution [-]')
ax_d.set_title(f'EOL decomposition (product = total)')
ax_d.grid(alpha=0.3); ax_d.legend(loc='lower left', fontsize=8)

fig.tight_layout()
fig.savefig(os.path.join(THIS_DIR, 'irradiation_scaling_factor.png'), dpi=140)
print("Plot written: irradiation_scaling_factor.png")

# ---------------------------------------------------------------------------
# Plot 4: GRSIS state at end of life (peak axial layer) — the raw model
# outputs behind the scaling-factor decomposition.
#   (a) swelling components: swclose is pinned at the interconnection
#       threshold s_th = 0.1 across the whole radius; swopen is small in the
#       hot-pressed inner core and spikes in the 2.9-3.3 mm band (too cool
#       for hot pressing, too hot for the alpha-phase sodium gate); swsol
#       (solid FP, 1.5%/at%) is radially flat and NOT part of the porosity
#       fed to the conductivity model.
#   (b) porosity split seen by the k model (legacy Baseir.for convention)
#       plus the sodium infiltration fraction f_inf.
#   (c) local Zr weight fraction from Zr redistribution — the driver of the
#       inner-core conductivity trend.
# ---------------------------------------------------------------------------
sol_eol = np.array([0.015 * bup_atpct[i_eol]] * NF)   # swsol (uniform, 1.5%/at%)
with h5py.File(OUTPUT_H5, "r") as f:
    sol_eol = f["grsis/swsol"][i_eol].reshape(NZ, NF)[j_peak]
    ptot_eol = f["thermal/poros_tot"][i_eol].reshape(NZ, NF)[j_peak]
    pgas_eol = f["thermal/poros_gas"][i_eol].reshape(NZ, NF)[j_peak]
    psod_eol = f["thermal/psod"][i_eol].reshape(NZ, NF)[j_peak]

fig, (axa, axb, axc) = plt.subplots(1, 3, figsize=(15, 4.6))
for a in (axa, axb, axc):
    a.axvline(0.6 * RFO * 1e3, color='k', ls=':', lw=1.0, alpha=0.7)
    a.axvline(r_phase_mm, color='g', ls='--', lw=1.0, alpha=0.7)
    a.set_xlabel('radius [mm]'); a.grid(alpha=0.3)

axa.step(r_mm, sc_eol, where='mid', label='$sw_{close}$ (closed bubbles)')
axa.step(r_mm, so_eol, where='mid', label='$sw_{open}$ (open bubbles)')
axa.step(r_mm, sol_eol, where='mid', label='$sw_{sol}$ (solid FP)')
axa.axhline(0.1, color='gray', lw=0.6, ls='--')
axa.text(0.05, 0.104, 'interconnection threshold $s_{th}$', fontsize=8, color='gray')
axa.set_ylabel('swelling fraction [-]')
axa.set_title('(a) GRSIS swelling components, EOL')
axa.legend(fontsize=8)

axb.step(r_mm, ptot_eol, where='mid', label='$p_{tot} = sw_{gas}$')
axb.step(r_mm, pgas_eol, where='mid', label='$p^{clos}_{gas} = sw_{close}$')
axb.step(r_mm, ptot_eol - pgas_eol, where='mid',
         label='open porosity $= sw_{open}$')
axb.step(r_mm, psod_eol, where='mid', ls='--',
         label='$f_{inf}$ (Na infiltration fraction)')
axb.set_ylabel('porosity / fraction [-]')
axb.set_title('(b) Porosity split fed to k model, EOL')
axb.legend(fontsize=8)

axc.step(r_mm, zr_eol, where='mid', color='C1')
axc.axhline(0.10, color='gray', lw=0.6, ls='--')
axc.text(0.05, 0.101, 'as-fabricated $x_{Zr}$', fontsize=8, color='gray')
axc.set_ylabel('Zr weight fraction [-]')
axc.set_title('(c) Zr redistribution, EOL')

fig.tight_layout()
fig.savefig(os.path.join(THIS_DIR, 'grsis_eol_profiles.png'), dpi=140)
print("Plot written: grsis_eol_profiles.png")

# Console summary used by legacy_fred_m_na/compare.py
print("\nSummary (peak axial layer, centerline):")
for t_target, lbl in [(1.0, '1 s'), (1000 * D2S, '1000 d'), (2176 * D2S, '2176 d')]:
    i = int(np.argmin(np.abs(times_s - t_target)))
    k_prof = k_fuel[i, j_peak*NF : (j_peak+1)*NF]
    print(f"  t={lbl:8s}  T_cl={T_fi_K[i, j_peak]:8.2f} K  "
          f"k(center)={k_prof[0]:.3f} W/m-K  k(edge)={k_prof[-1]:.3f} W/m-K  "
          f"bup={bup_atpct[i]:.2f} at%")