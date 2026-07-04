"""
FRED-M-Na: Timpano/SAS4A base-irradiation benchmark
U-Pu-Zr metallic fuel (13.25 wt% Pu, 10 wt% Zr) / HT-9 cladding
Sodium-cooled fast reactor, based on ESFR-SIMPLE geometry

Input parameters from ref_input.txt.
Compared against: legacy FRED-M and SAS4A at t=0, 45, 135, 365, 730, 1088, 2176 days.

Outputs:
  - fred_m_na_base_irradiation.h5 : full simulation results (HDF5, written by C++ at each step)
  - comparison_verification.png   : 6-panel figure comparing temperatures + extra variables
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
# Power table (from ref_input.txt, card 300000)
# Format: [time_s, qqv_z1, qqv_z2, ..., qqv_z24]  [W/m3]
# ---------------------------------------------------------------------------
POWER_ROWS = [
    [0.0000E+00, 7.1054E+08, 7.3022E+08, 7.6210E+08, 8.0707E+08, 8.5014E+08, 8.9091E+08, 9.2669E+08, 9.5881E+08, 9.8458E+08, 1.0031E+09, 1.0172E+09, 1.0220E+09, 1.0212E+09, 1.0148E+09, 9.9911E+08, 9.7890E+08, 9.5106E+08, 9.1658E+08, 8.7825E+08, 8.3429E+08, 7.8734E+08, 7.3771E+08, 6.9631E+08, 6.6843E+08],
    [3.9096E+06, 7.0953E+08, 7.2918E+08, 7.6102E+08, 8.0593E+08, 8.4894E+08, 8.8965E+08, 9.2537E+08, 9.5745E+08, 9.8319E+08, 1.0017E+09, 1.0157E+09, 1.0205E+09, 1.0198E+09, 1.0133E+09, 9.9769E+08, 9.7751E+08, 9.4971E+08, 9.1528E+08, 8.7701E+08, 8.3311E+08, 7.8622E+08, 7.3667E+08, 6.9532E+08, 6.6748E+08],
    [7.8192E+06, 7.0852E+08, 7.2815E+08, 7.5994E+08, 8.0478E+08, 8.4773E+08, 8.8839E+08, 9.2406E+08, 9.5609E+08, 9.8179E+08, 1.0003E+09, 1.0143E+09, 1.0191E+09, 1.0183E+09, 1.0119E+09, 9.9628E+08, 9.7612E+08, 9.4836E+08, 9.1398E+08, 8.7576E+08, 8.3192E+08, 7.8511E+08, 7.3562E+08, 6.9434E+08, 6.6653E+08],
    [1.1729E+07, 7.0752E+08, 7.2711E+08, 7.5886E+08, 8.0364E+08, 8.4652E+08, 8.8712E+08, 9.2274E+08, 9.5473E+08, 9.8039E+08, 9.9883E+08, 1.0128E+09, 1.0176E+09, 1.0169E+09, 1.0105E+09, 9.9486E+08, 9.7473E+08, 9.4701E+08, 9.1268E+08, 8.7451E+08, 8.3074E+08, 7.8399E+08, 7.3458E+08, 6.9335E+08, 6.6558E+08],
    [1.5638E+07, 7.0651E+08, 7.2608E+08, 7.5777E+08, 8.0249E+08, 8.4532E+08, 8.8586E+08, 9.2143E+08, 9.5337E+08, 9.7900E+08, 9.9741E+08, 1.0114E+09, 1.0162E+09, 1.0154E+09, 1.0090E+09, 9.9344E+08, 9.7334E+08, 9.4566E+08, 9.1138E+08, 8.7327E+08, 8.2956E+08, 7.8287E+08, 7.3353E+08, 6.9236E+08, 6.6464E+08],
    [1.9548E+07, 7.0550E+08, 7.2504E+08, 7.5669E+08, 8.0135E+08, 8.4411E+08, 8.8460E+08, 9.2012E+08, 9.5201E+08, 9.7760E+08, 9.9598E+08, 1.0100E+09, 1.0147E+09, 1.0140E+09, 1.0076E+09, 9.9202E+08, 9.7196E+08, 9.4431E+08, 9.1008E+08, 8.7202E+08, 8.2837E+08, 7.8176E+08, 7.3248E+08, 6.9137E+08, 6.6369E+08],
    [2.3458E+07, 7.0449E+08, 7.2401E+08, 7.5561E+08, 8.0020E+08, 8.4291E+08, 8.8333E+08, 9.1880E+08, 9.5065E+08, 9.7620E+08, 9.9456E+08, 1.0085E+09, 1.0133E+09, 1.0125E+09, 1.0061E+09, 9.9061E+08, 9.7057E+08, 9.4296E+08, 9.0878E+08, 8.7078E+08, 8.2719E+08, 7.8064E+08, 7.3144E+08, 6.9039E+08, 6.6274E+08],
    [2.7367E+07, 7.0348E+08, 7.2297E+08, 7.5453E+08, 7.9906E+08, 8.4170E+08, 8.8207E+08, 9.1749E+08, 9.4929E+08, 9.7481E+08, 9.9314E+08, 1.0071E+09, 1.0118E+09, 1.0111E+09, 1.0047E+09, 9.8919E+08, 9.6918E+08, 9.4161E+08, 9.0748E+08, 8.6953E+08, 8.2601E+08, 7.7952E+08, 7.3039E+08, 6.8940E+08, 6.6179E+08],
    [3.1277E+07, 7.0248E+08, 7.2194E+08, 7.5345E+08, 7.9791E+08, 8.4049E+08, 8.8081E+08, 9.1617E+08, 9.4793E+08, 9.7341E+08, 9.9172E+08, 1.0056E+09, 1.0104E+09, 1.0096E+09, 1.0033E+09, 9.8777E+08, 9.6779E+08, 9.4027E+08, 9.0618E+08, 8.6829E+08, 8.2482E+08, 7.7841E+08, 7.2934E+08, 6.8841E+08, 6.6084E+08],
    [3.1363E+07, 7.1296E+08, 7.3271E+08, 7.6470E+08, 8.0982E+08, 8.5304E+08, 8.9395E+08, 9.2985E+08, 9.6208E+08, 9.8794E+08, 1.0065E+09, 1.0206E+09, 1.0255E+09, 1.0247E+09, 1.0182E+09, 1.0025E+09, 9.8223E+08, 9.5430E+08, 9.1970E+08, 8.8124E+08, 8.3713E+08, 7.9002E+08, 7.4023E+08, 6.9868E+08, 6.7071E+08],
    [4.7002E+07, 7.0772E+08, 7.2733E+08, 7.5908E+08, 8.0388E+08, 8.4677E+08, 8.8739E+08, 9.2302E+08, 9.5502E+08, 9.8068E+08, 9.9913E+08, 1.0131E+09, 1.0179E+09, 1.0172E+09, 1.0108E+09, 9.9515E+08, 9.7502E+08, 9.4729E+08, 9.1295E+08, 8.7477E+08, 8.3099E+08, 7.8422E+08, 7.3479E+08, 6.9355E+08, 6.6578E+08],
    [6.2640E+07, 7.0249E+08, 7.2195E+08, 7.5347E+08, 7.9793E+08, 8.4051E+08, 8.8082E+08, 9.1619E+08, 9.4795E+08, 9.7343E+08, 9.9173E+08, 1.0056E+09, 1.0104E+09, 1.0096E+09, 1.0033E+09, 9.8779E+08, 9.6781E+08, 9.4028E+08, 9.0619E+08, 8.6830E+08, 8.2484E+08, 7.7842E+08, 7.2936E+08, 6.8842E+08, 6.6086E+08],
    [6.2726E+07, 7.1520E+08, 7.3501E+08, 7.6710E+08, 8.1237E+08, 8.5572E+08, 8.9676E+08, 9.3277E+08, 9.6511E+08, 9.9105E+08, 1.0097E+09, 1.0238E+09, 1.0287E+09, 1.0279E+09, 1.0214E+09, 1.0057E+09, 9.8532E+08, 9.5730E+08, 9.2259E+08, 8.8402E+08, 8.3977E+08, 7.9251E+08, 7.4256E+08, 7.0088E+08, 6.7282E+08],
    [7.8365E+07, 7.1260E+08, 7.3234E+08, 7.6431E+08, 8.0941E+08, 8.5261E+08, 8.9350E+08, 9.2938E+08, 9.6160E+08, 9.8744E+08, 1.0060E+09, 1.0201E+09, 1.0250E+09, 1.0242E+09, 1.0177E+09, 1.0020E+09, 9.8174E+08, 9.5382E+08, 9.1924E+08, 8.8080E+08, 8.3671E+08, 7.8963E+08, 7.3986E+08, 6.9833E+08, 6.7037E+08],
    [9.4003E+07, 7.1000E+08, 7.2967E+08, 7.6152E+08, 8.0646E+08, 8.4950E+08, 8.9024E+08, 9.2598E+08, 9.5809E+08, 9.8384E+08, 1.0023E+09, 1.0164E+09, 1.0212E+09, 1.0204E+09, 1.0140E+09, 9.9835E+08, 9.7816E+08, 9.5034E+08, 9.1588E+08, 8.7758E+08, 8.3366E+08, 7.8674E+08, 7.3715E+08, 6.9578E+08, 6.6792E+08],
    [9.4090E+07, 7.1800E+08, 7.3788E+08, 7.7010E+08, 8.1554E+08, 8.5906E+08, 9.0026E+08, 9.3641E+08, 9.6888E+08, 9.9492E+08, 1.0136E+09, 1.0278E+09, 1.0327E+09, 1.0319E+09, 1.0254E+09, 1.0096E+09, 9.8917E+08, 9.6104E+08, 9.2619E+08, 8.8747E+08, 8.4305E+08, 7.9560E+08, 7.4546E+08, 7.0362E+08, 6.7544E+08],
    [1.0973E+08, 7.1206E+08, 7.3179E+08, 7.6373E+08, 8.0880E+08, 8.5197E+08, 8.9283E+08, 9.2868E+08, 9.6087E+08, 9.8670E+08, 1.0052E+09, 1.0194E+09, 1.0242E+09, 1.0234E+09, 1.0170E+09, 1.0013E+09, 9.8100E+08, 9.5310E+08, 9.1854E+08, 8.8014E+08, 8.3608E+08, 7.8903E+08, 7.3930E+08, 6.9780E+08, 6.6986E+08],
    [1.2537E+08, 7.0613E+08, 7.2569E+08, 7.5737E+08, 8.0207E+08, 8.4487E+08, 8.8539E+08, 9.2094E+08, 9.5286E+08, 9.7847E+08, 9.9687E+08, 1.0109E+09, 1.0157E+09, 1.0149E+09, 1.0085E+09, 9.9291E+08, 9.7283E+08, 9.4516E+08, 9.1089E+08, 8.7280E+08, 8.2911E+08, 7.8246E+08, 7.3314E+08, 6.9199E+08, 6.6428E+08],
    [1.2545E+08, 7.1523E+08, 7.3504E+08, 7.6713E+08, 8.1240E+08, 8.5576E+08, 8.9680E+08, 9.3281E+08, 9.6515E+08, 9.9108E+08, 1.0097E+09, 1.0239E+09, 1.0287E+09, 1.0280E+09, 1.0215E+09, 1.0057E+09, 9.8536E+08, 9.5734E+08, 9.2263E+08, 8.8405E+08, 8.3980E+08, 7.9254E+08, 7.4259E+08, 7.0091E+08, 6.7284E+08],
    [1.4109E+08, 7.0714E+08, 7.2673E+08, 7.5846E+08, 8.0322E+08, 8.4608E+08, 8.8666E+08, 9.2226E+08, 9.5423E+08, 9.7988E+08, 9.9831E+08, 1.0123E+09, 1.0171E+09, 1.0163E+09, 1.0099E+09, 9.9434E+08, 9.7422E+08, 9.4651E+08, 9.1220E+08, 8.7406E+08, 8.3031E+08, 7.8358E+08, 7.3419E+08, 6.9299E+08, 6.6523E+08],
    [1.5673E+08, 6.9906E+08, 7.1842E+08, 7.4979E+08, 7.9403E+08, 8.3641E+08, 8.7652E+08, 9.1172E+08, 9.4332E+08, 9.6868E+08, 9.8689E+08, 1.0007E+09, 1.0055E+09, 1.0047E+09, 9.9839E+08, 9.8297E+08, 9.6308E+08, 9.3569E+08, 9.0177E+08, 8.6406E+08, 8.2081E+08, 7.7462E+08, 7.2580E+08, 6.8506E+08, 6.5763E+08],
    [1.5682E+08, 7.1264E+08, 7.3238E+08, 7.6436E+08, 8.0946E+08, 8.5266E+08, 8.9355E+08, 9.2943E+08, 9.6165E+08, 9.8750E+08, 1.0061E+09, 1.0202E+09, 1.0250E+09, 1.0242E+09, 1.0178E+09, 1.0021E+09, 9.8180E+08, 9.5387E+08, 9.1929E+08, 8.8085E+08, 8.3676E+08, 7.8967E+08, 7.3990E+08, 6.9837E+08, 6.7041E+08],
    [1.7245E+08, 7.0695E+08, 7.2653E+08, 7.5825E+08, 8.0300E+08, 8.4585E+08, 8.8642E+08, 9.2201E+08, 9.5397E+08, 9.7961E+08, 9.9803E+08, 1.0120E+09, 1.0168E+09, 1.0161E+09, 1.0097E+09, 9.9407E+08, 9.7396E+08, 9.4626E+08, 9.1195E+08, 8.7382E+08, 8.3008E+08, 7.8337E+08, 7.3399E+08, 6.9280E+08, 6.6505E+08],
    [1.8809E+08, 7.0126E+08, 7.2068E+08, 7.5215E+08, 7.9653E+08, 8.3904E+08, 8.7928E+08, 9.1459E+08, 9.4629E+08, 9.7173E+08, 9.9000E+08, 1.0039E+09, 1.0086E+09, 1.0079E+09, 1.0015E+09, 9.8606E+08, 9.6612E+08, 9.3864E+08, 9.0461E+08, 8.6678E+08, 8.2340E+08, 7.7706E+08, 7.2808E+08, 6.8722E+08, 6.5970E+08],
]
POWER_ROWS = np.array(POWER_ROWS)
power_times = POWER_ROWS[:, 0].tolist()
power_qqv   = POWER_ROWS[:, 1:].T.tolist()
# Constant extrapolation past end
power_times.append(2.0e8)
for j in range(24):
    power_qqv[j].append(power_qqv[j][-1])

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
solver.set_conductivity_model(fred.ConductivityModel.EsfrSimple)
solver.set_power_density_history_per_layer(power_times, power_qqv)
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
# no SUNDIALS IDA — see FredMNaSolver.hpp) needs an explicit internal step
# size: DTOUT below (90 d) is only the output-writing cadence, and taking
# that as one giant physics step would under-resolve the threshold-driven
# fanis "soft" contact transition (~day 100-130) and other fast-varying
# quantities (gpres, FGR) — see scripts/session_state.md's step-size
# comparison. 1 day balances accuracy against wall-clock time; unlike the
# old SUNDIALS IDA-based scheme this replaced, there is no adaptive-step
# collapse or reinit-failure risk here regardless of step size choice.
solver.set_step_size(86400.0)

# ---------------------------------------------------------------------------
# Run: full depletion 0 → 2176 days, with HDF5 output (C++ streaming)
# ---------------------------------------------------------------------------
D2S = 86400.0
TEND = 2176 * D2S
DTOUT = 90  * D2S   # output every 90 days (~24 steps over full run)

THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
OUTPUT_H5  = os.path.join(THIS_DIR, "fred_m_na_base_irradiation.h5")
HOTSTART_H5 = os.path.join(THIS_DIR, "fred_m_na_hotstart.h5")

print("=" * 70)
print("FRED-M-Na: Timpano SAS4A Base Irradiation Benchmark")
print("=" * 70)
print(f"  Fuel: U-Pu-Zr  Pu={fuel.pu_content():.4f}  Zr={fuel.zr_content():.4f}")
print(f"  rfi={RFI*1e3:.3f} mm  rfo={RFO*1e3:.3f} mm  rci={RCI*1e3:.3f} mm  rco={RCO*1e3:.3f} mm")
print(f"  nf={NF}  nc={NC}  nz={NZ}  dz={DZ0*1e3:.1f} mm  active height={NZ*DZ0:.3f} m")
print(f"  Conductivity: EsfrSimple (f=3 ESFR-SIMPLE sigmoid)  |  HTC: Subbotin")
print(f"  tend={TEND/D2S:.0f} d  dtout={DTOUT/D2S:.0f} d  → HDF5: {os.path.basename(OUTPUT_H5)}")
print()

# Single 1-day step from cold start → thermal equilibrium with negligible burnup.
# Used for the "beginning of life" temperature comparison panel instead of the
# literal t=0 (cold, 633 K flat) or t=90 d (first main output, carries some burnup).
solver.set_output_file(HOTSTART_H5)
solver.run(1 * D2S, 1 * D2S)
print("Hot-start capture complete (t=1 d).")

# Full irradiation run (resets to t=0 internally, independent of the above).
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
    gpres_mna = f["burnup/gpres"][:]          # [MPa]
    xwast_mna = f["burnup/xwast"][:]          # [m] max across layers
    fggen_mna = f["burnup/fggen"][:]          # [mol] total generated
    fgrel_mna = f["burnup/fgrel"][:]          # [mol] total released

times_d = times_s / D2S
stride  = NF + NC

# Fuel inner (centerline) temperature per axial layer, shape (n_steps, NZ)
T_fi_K = T_flat[:, 0::stride]
T_fi_C = T_fi_K - 273.15

# Hot-start snapshot (t=1 d): thermal equilibrium, negligible burnup.
with h5py.File(HOTSTART_H5, "r") as f:
    T_flat_hot = f["thermal/T"][:]   # shape (1, nz*(nf+nc))
T_fi_C_hot = T_flat_hot[:, 0::stride] - 273.15   # shape (1, NZ)

# FGR fraction
fgr_frac_mna = np.zeros_like(fggen_mna)
_mask = fggen_mna > 0
fgr_frac_mna[_mask] = fgrel_mna[_mask] / fggen_mna[_mask]

# Print key result table
print(f"\n{'t [d]':>8} {'gpres [MPa]':>12} {'xwast_max [µm]':>16} {'FGR [%]':>10} {'T_fi_peak [°C]':>16}")
print("-" * 70)
for i, td in enumerate(times_d):
    T_peak_C = T_fi_C[i].max()
    print(f"{td:8.1f} {gpres_mna[i]:12.4f} {xwast_mna[i]*1e6:16.1f} "
          f"{fgr_frac_mna[i]*100:10.2f} {T_peak_C:16.1f}")

# ---------------------------------------------------------------------------
# Temperature comparison vs reference data
# ---------------------------------------------------------------------------
try:
    from reference_thermal_data import data as ref_data
    HAS_REF = True
except ImportError:
    HAS_REF = False
    print("WARNING: reference_thermal_data.py not found — skipping T comparison")

try:
    from extra_variables_fred_legacy import snapshots as legacy_snap
    HAS_LEGACY = True
except ImportError:
    HAS_LEGACY = False
    print("WARNING: extra_variables_fred_legacy.py not found — skipping extra variable comparison")

z_centers = np.array([0.01562 + j * 0.03125 for j in range(NZ)])

# ---------------------------------------------------------------------------
# Figure: 6-panel comparison (3 temperature + 3 extra variables)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle(
    "FRED-M-Na vs Legacy FRED-M vs SAS4A — Timpano SAS4A Benchmark\n"
    f"U-Pu-Zr (13.25%Pu / 10%Zr) HT-9 clad | Na coolant | EsfrSimple + Subbotin",
    fontsize=11
)

# ---- Row 1: Fuel inner temperature vs z at 3 reference times ----
TARGET_DAYS = [0, 1088, 2176]
for col, tgt_day in enumerate(TARGET_DAYS):
    ax = axes[0, col]
    if tgt_day == 0:
        # FRED-M-Na starts cold; the dedicated hot-start run captures t=1 d
        # (thermal equilibrium, negligible burnup) as the BOL reference.
        ax.plot(z_centers, T_fi_C_hot[0], 'b-o', ms=3, lw=1.5,
                label='FRED-M-Na (t=1 d, hot steady state)')
        t_mna = 1.0
    else:
        idx_mna = np.argmin(np.abs(times_d - tgt_day))
        t_mna = times_d[idx_mna]
        ax.plot(z_centers, T_fi_C[idx_mna], 'b-o', ms=3, lw=1.5,
                label=f'FRED-M-Na (t={t_mna:.0f} d)')

    if HAS_REF and tgt_day in ref_data:
        fred_r = ref_data[tgt_day].get('FRED', {})
        T_fi_fred = fred_r.get('Tfuel_inner(C)', [])
        z_fred    = fred_r.get('z (m)', [])
        if T_fi_fred and len(T_fi_fred) == NZ:
            ax.plot(z_fred, T_fi_fred, 'r--s', ms=3, lw=1.5, label='Legacy FRED-M')

        sas_r = ref_data[tgt_day].get('SAS4A', {})
        T_fi_sas = sas_r.get('Tfuel_inner(C)', [])
        z_sas    = sas_r.get('z (m)', [])
        if T_fi_sas and len(T_fi_sas) == NZ:
            order = np.argsort(z_sas)
            ax.plot(np.array(z_sas)[order], np.array(T_fi_sas)[order],
                    'g:^', ms=3, lw=1.5, label='SAS4A')

    ax.set_xlabel('Axial position z [m]', fontsize=9)
    ax.set_ylabel('T fuel inner [°C]', fontsize=9)
    ax.set_title(f't = {tgt_day} days', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Print numerical comparison
    if HAS_REF and tgt_day in ref_data:
        fred_r = ref_data[tgt_day].get('FRED', {})
        T_fi_fred = fred_r.get('Tfuel_inner(C)', [])
        if T_fi_fred and len(T_fi_fred) == NZ:
            T_fi_fred_K = np.array(T_fi_fred) + 273.15
            T_mna_K = (T_fi_C_hot[0] + 273.15) if tgt_day == 0 else T_fi_K[idx_mna]
            err = T_mna_K - T_fi_fred_K
            print(f"\nt={tgt_day} d (using FRED-M-Na t={t_mna:.0f} d):")
            print(f"  vs FRED: dT_fi mean={np.mean(err):+.1f} K  max={np.max(np.abs(err)):.1f} K")

# ---- Row 2: Extra variable comparisons vs time ----
if HAS_LEGACY:
    leg_keys = sorted(legacy_snap)
    leg_t         = np.array([legacy_snap[k]['time_d'] for k in leg_keys])
    leg_gpres     = np.array([legacy_snap[k]['gpres_MPa'] for k in leg_keys])
    leg_xwast_max = np.array([legacy_snap[k]['xwast_um'].max() for k in leg_keys])
    leg_fgr_max   = np.array([legacy_snap[k]['fgr_frac'].max() for k in leg_keys])

    # Plenum pressure
    ax = axes[1, 0]
    ax.plot(leg_t, leg_gpres, 'r-o', ms=5, lw=1.5, label='Legacy FRED-M')
    ax.plot(times_d, gpres_mna, 'b-^', ms=5, lw=1.5, label='FRED-M-Na')
    ax.set_xlabel('Time [days]', fontsize=9)
    ax.set_ylabel('Plenum gas pressure [MPa]', fontsize=9)
    ax.set_title('Gas Pressure', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Max cladding wastage
    ax = axes[1, 1]
    ax.plot(leg_t, leg_xwast_max, 'r-o', ms=5, lw=1.5, label='Legacy FRED-M')
    ax.plot(times_d, xwast_mna * 1e6, 'b-^', ms=5, lw=1.5, label='FRED-M-Na')
    ax.set_xlabel('Time [days]', fontsize=9)
    ax.set_ylabel('Max cladding wastage [µm]', fontsize=9)
    ax.set_title('Lanthanide Cladding Wastage', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # FGR fraction (max axial node for legacy; fgrel/fggen for FRED-M-Na)
    ax = axes[1, 2]
    ax.plot(leg_t, leg_fgr_max * 100, 'r-o', ms=5, lw=1.5,
            label='Legacy FRED-M\n(rod-total fgrel/fggen)')
    ax.plot(times_d, fgr_frac_mna * 100, 'b-^', ms=5, lw=1.5,
            label='FRED-M-Na\n(total fgrel/fggen)')
    ax.set_xlabel('Time [days]', fontsize=9)
    ax.set_ylabel('FGR fraction [%]', fontsize=9)
    ax.set_title('Fission Gas Release', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Print comparison at final time
    print("\nExtra variable comparison at end of run:")
    print(f"  gpres:   Legacy={leg_gpres[-1]:.3f} MPa  FRED-M-Na={gpres_mna[-1]:.3f} MPa")
    print(f"  xwast:   Legacy={leg_xwast_max[-1]:.1f} µm   FRED-M-Na={xwast_mna[-1]*1e6:.1f} µm")
    print(f"  fgr:     Legacy={leg_fgr_max[-1]*100:.2f}%   FRED-M-Na={fgr_frac_mna[-1]*100:.2f}%")
else:
    for col in range(3):
        axes[1, col].text(0.5, 0.5, 'extra_variables_fred_legacy.py not found',
                          ha='center', va='center', transform=axes[1, col].transAxes)

plt.tight_layout(rect=[0, 0, 1, 0.95])
OUT_PNG = os.path.join(THIS_DIR, "comparison_verification.png")
fig.savefig(OUT_PNG, dpi=150, bbox_inches='tight')
print(f"\nSaved: {os.path.basename(OUT_PNG)}")
plt.close(fig)
print("Done.")
