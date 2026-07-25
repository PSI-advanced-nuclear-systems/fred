#!/usr/bin/env python3
"""
Cladding wastage — LaTracking model (MFUEL App. 9.3 lanthanide diffusion).

Per-fuel-node radial diffusion of the lanthanide number density u(r,t) with
a uniform fission source; while the gap is in soft/clos contact the
fuel-cladding interface acts as a Dirichlet sink and the time-integrated
arriving flux, divided by the HT-9 uptake capacity u_sol, gives the wastage
depth.  With an open gap the interface is zero-flux: lanthanides pile up in
the fuel and no wastage grows — the threshold-demo plot shows this.

Outputs: fred_m_na_la.h5 (incl. burnup/xwast_layer and per-node burnup/c_la),
         wastage_threshold_demo.png, la_radial_profiles.png
"""

import sys, os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'build'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import _fred_m_na as fred

from common import (make_solver, plot_threshold_demo, pick_demo_layer,
                    read_wastage, load_benchmark_power,
                    TEND, DTOUT, D2S, NF, NZ, RFO)

THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
OUTPUT_H5 = os.path.join(THIS_DIR, "fred_m_na_la.h5")

solver = make_solver(fred)
solver.set_cladding_wastage_model(fred.CladWastageModel.LaTracking)
solver.set_output_file(OUTPUT_H5)
solver.run(TEND, DTOUT)
print("Run complete.")

ok = plot_threshold_demo(OUTPUT_H5,
                         os.path.join(THIS_DIR, 'wastage_threshold_demo.png'),
                         'LaTracking: wastage gated on soft/clos contact')

# ---------------------------------------------------------------------------
# Radial La profiles at the demo layer: BOL, just before first soft contact
# (open gap: zero-flux boundary, La piles up at the rim), shortly after soft
# contact (sink active: rim clamped to zero, near-edge gradient develops),
# after clos, and EOL (quasi-steady source/sink balance).
# ---------------------------------------------------------------------------
t, bu, xw, con = read_wastage(OUTPUT_H5)
with h5py.File(OUTPUT_H5, 'r') as f:
    c_la = f["burnup/c_la"][:]          # (n, NZ*NF). n = # of time steps. see reshape c_la[i].reshape(NZ, NF)

j, _ = pick_demo_layer(t, con)
i_soft = np.nonzero(con[:, j] >= 1.0)[0][0]
i_clos_all = np.nonzero(con[:, j] >= 2.0)[0]
i_clos = i_clos_all[0] if i_clos_all.size else None

picks = [(0,                    't = 1 s (BOL)'),
         (max(i_soft - 1, 0),   f'just before soft ({t[i_soft-1]/D2S:.0f} d, open gap)'),
         (min(i_soft + 30, len(t)-1),
                                f'soft + 30 d ({t[min(i_soft+30, len(t)-1)]/D2S:.0f} d, sink on)')]
if i_clos is not None:
    picks.append((min(i_clos + 100, len(t)-1),
                  f'clos + 100 d ({t[min(i_clos+100, len(t)-1)]/D2S:.0f} d)'))
picks.append((len(t) - 1, f't = {t[-1]/D2S:.0f} d (EOL)'))

r_mm = np.linspace(0.0, RFO, NF) * 1e3
fig, ax = plt.subplots(figsize=(8.5, 5))
for i, lbl in picks:
    prof = c_la[i].reshape(NZ, NF)[j]
    ax.plot(r_mm, prof, 'o-', ms=3.5, label=lbl)
ax.set_xlabel('radius [mm]')
ax.set_ylabel('lanthanide number density $u$ [#/m³]')
ax.set_title(f'LaTracking: radial La profiles (layer {j+1}/{NZ})')
ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout()
out = os.path.join(THIS_DIR, 'la_radial_profiles.png')
fig.savefig(out, dpi=140)
print(f"Plot written: {out}")

# ---------------------------------------------------------------------------
# Lanthanide conservation check at the demo layer:
#   generated (Y_La * integral q'''(t) dt * V) = retained in fuel
#   + absorbed by cladding (absorbed = xwast * u_sol * interface area)
# ---------------------------------------------------------------------------
Y_LA, U_SOL = 15.29e9, 6.4e27
dr = RFO / (NF - 1)
V = np.zeros(NF)                      # control volumes per unit length / 2*pi
V[0] = 0.5 * (0.5*dr)**2
for i in range(1, NF-1): V[i] = (i*dr) * dr
V[-1] = ((NF-1)*dr - 0.25*dr) * (0.5*dr)

# Time-integrated power density of the demo layer from the benchmark table
p_times, p_qqv = load_benchmark_power()
tt = np.asarray(p_times); qq = np.asarray(p_qqv[j])
tt_clip = np.clip(tt, 0.0, t[-1])
E_dep = np.trapezoid(qq, tt_clip)                 # [J/m3] deposited to EOL

gen      = Y_LA * E_dep * V.sum()                 # per unit length / 2*pi
retained = float((c_la[-1].reshape(NZ, NF)[j] * V).sum())
absorbed = xw[-1, j] * U_SOL * RFO                # per unit length / 2*pi
balance  = (retained + absorbed) / gen
print(f"La conservation (layer {j+1}): generated={gen:.4e}  retained={retained:.4e}"
      f"  absorbed={absorbed:.4e}  (retained+absorbed)/generated = {balance:.4f}")
# The solver's subgrid solve (la_refine=8 by default) is exactly conservative;
# this check reconstructs the retained inventory from the NF coarse samples in
# the H5, which cannot resolve the sub-cell interface boundary layer, so a
# few-percent quadrature error is expected here (it is NOT a solver leak).
TOL = 0.05
ok = ok and abs(balance - 1.0) < TOL
print("Conservation check (coarse-sampled, tol 5%):",
      "PASS" if abs(balance - 1.0) < TOL else "FAIL")

sys.exit(0 if ok else 1)
