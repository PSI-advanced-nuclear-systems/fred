#!/usr/bin/env python3
"""
Compare FRED platform (../run.py, DetailedNaSodium) against legacy FRED-M
(fred.x with FUEL_THK 1) for the sodium-infiltration case:

  - peak fuel centerline temperature (axial layer 24, the hottest under
    axially uniform power) at ~1 s (BOL), 1000 d, and 2176 d
  - fuel thermal conductivity radial profile at the same times

Inputs
  ../fred_m_na_base_irradiation.h5 : platform results (thermal/T, thermal/k_fuel)
  fred.dat                         : legacy per-output scalars at OUTPUT_AXIAL_LAYER 24
  kfuel.dat                        : legacy kfuel_m(i,j) all layers, per output event

Grid note: legacy kfuel_m lives on the nf-1 radial interval midpoints
(T = 0.5*(tem_i + tem_i+1)); the platform saves nf node values.  For the
comparison the platform profile is averaged to midpoints: 0.5*(k_i + k_i+1).

Legacy's first output is at t = nout*dtout = 1e5 s (~1.16 d), so the "1 s"
column compares the platform's t = 1 s hot-start state against legacy's
first available output — both are at thermal equilibrium with negligible
burnup (< 0.01 at%).
"""

import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
H5   = os.path.join(THIS_DIR, '..', 'fred_m_na_base_irradiation.h5')
NF, NC, NZ = 11, 5, 24
J = 23          # axial layer 24 (0-based), OUTPUT_AXIAL_LAYER 24 in the deck
D2S = 86400.0
RFO = 3.863e-3

# ── Platform results ─────────────────────────────────────────────────────────
with h5py.File(H5, 'r') as f:
    p_t  = f['time'][:]
    p_T  = f['thermal/T'][:]            # (n, NZ*(NF+NC))
    p_k  = f['thermal/k_fuel'][:]       # (n, NZ*NF)
    p_bu = f['burnup/bup_atpct'][:]

stride  = NF + NC
p_Tcl   = p_T[:, J*stride]                       # layer-24 centerline [K]
p_kmid  = 0.5 * (p_k[:, J*NF:(J+1)*NF-1] + p_k[:, J*NF+1:(J+1)*NF])  # midpoints

# ── Legacy results ───────────────────────────────────────────────────────────
with open(os.path.join(THIS_DIR, 'fred.dat')) as f:
    header = f.readline().split()
cols = {name: i for i, name in enumerate(header)}
raw  = np.genfromtxt(os.path.join(THIS_DIR, 'fred.dat'), skip_header=1,
                     usecols=range(len(header) - 1))   # drop gap_state string
l_t   = raw[:, cols['time(s)']]
l_Tcl = raw[:, cols['temfi(c)']] + 273.15              # layer-24 centerline [K]
# NOTE: fred.dat's bupave(FIMA) column is legacy's *Yichao/Olander* burnup
# metric (205 MeV/fission, average molar mass with a (1-zr) correction).
# Legacy's PHYSICS (sinf/flamb calls in Baseir.for) uses bup_FIMA2, which is
# formula-identical to the platform's (200 MeV, HM atom density).  Rescale
# the reported column by the exact constant-composition ratio so the burnup
# rows compare like-for-like; at equal power the two codes' physics burnups
# are then identical by construction.
_pu, _zr = 0.1325, 0.10
_u  = 1.0 - _pu - _zr
_M  = _zr*0.091224 + _pu*0.244 + _u*0.238029       # avg molar mass [kg/mol]
_S  = _pu/0.244 + _u/0.238029                      # HM mol per kg fuel [mol/kg]
_bu_scale = (3.28451e-11 / 3.204354e-11) * (1.0 - _zr) / (_M * _S)
l_bu  = raw[:, cols['bupave(FIMA)']] * 100.0 * _bu_scale   # -> platform metric [at%]

kraw   = np.loadtxt(os.path.join(THIS_DIR, 'kfuel.dat'), comments='!')
sel    = kraw[:, 2].astype(int) == J + 1               # layer-24 rows
lk_t   = kraw[sel, 0]
lk_k   = kraw[sel, 3:3 + NF - 1]                       # interval midpoints

# ── Comparison table ─────────────────────────────────────────────────────────
r_mid_mm = 0.5 * (np.linspace(0, RFO, NF)[:-1] + np.linspace(0, RFO, NF)[1:]) * 1e3

print("=" * 78)
print("Platform (DetailedNaSodium) vs legacy FRED-M (FUEL_THK 1) — layer 24")
print("=" * 78)
rows = [(1.0, "1 s (BOL)"), (1088 * D2S, "1088 d (mid-life)"), (2176 * D2S, "2176 d")]
for t_target, lbl in rows:
    ip = int(np.argmin(np.abs(p_t - t_target)))
    il = int(np.argmin(np.abs(l_t - t_target)))
    ik = int(np.argmin(np.abs(lk_t - t_target)))
    kp, kl = p_kmid[ip], lk_k[ik]
    dk = 100.0 * np.max(np.abs(kp - kl) / kl)
    print(f"\n  {lbl}   (platform t={p_t[ip]:.3g} s, legacy t={l_t[il]:.3g} s)")
    print(f"    T_centerline : platform {p_Tcl[ip]:8.2f} K | legacy {l_Tcl[il]:8.2f} K "
          f"| diff {p_Tcl[ip]-l_Tcl[il]:+7.2f} K")
    print(f"    burnup [at%] : platform {p_bu[ip]:8.2f}   | legacy {l_bu[il]:8.2f}")
    print(f"    k centre     : platform {kp[0]:8.3f}   | legacy {kl[0]:8.3f}   W/m-K")
    print(f"    k edge       : platform {kp[-1]:8.3f}   | legacy {kl[-1]:8.3f}   W/m-K")
    print(f"    max |dk|/k over radius: {dk:.2f} %")

# ── Plots ────────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

ax1.plot(p_t / D2S, p_Tcl - 273.15, '-',  label='FRED platform (DetailedNaSodium)')
ax1.plot(l_t / D2S, l_Tcl - 273.15, '--', label='legacy FRED-M (FUEL_THK 1)')
ax1.set_xlabel('time [days]'); ax1.set_ylabel('centerline T, layer 24 [°C]')
ax1.grid(alpha=0.3); ax1.legend(); ax1.set_title('Peak fuel centerline temperature')

colors = ['C0', 'C1', 'C2']
for (t_target, lbl), c in zip(rows, colors):
    ip = int(np.argmin(np.abs(p_t - t_target)))
    ik = int(np.argmin(np.abs(lk_t - t_target)))
    ax2.plot(r_mid_mm, p_kmid[ip], 'o-',  color=c, label=f'platform, {lbl}')
    ax2.plot(r_mid_mm, lk_k[ik],   's--', color=c, label=f'legacy, {lbl}')
ax2.set_xlabel('radius [mm]'); ax2.set_ylabel('fuel k [W/m·K]')
ax2.grid(alpha=0.3); ax2.legend(fontsize=8)
ax2.set_title('Fuel conductivity, radial profile (layer 24)')

fig.tight_layout()
out = os.path.join(THIS_DIR, 'comparison_platform_vs_legacy.png')
fig.savefig(out, dpi=140)
print(f"\nPlot written: {out}")
