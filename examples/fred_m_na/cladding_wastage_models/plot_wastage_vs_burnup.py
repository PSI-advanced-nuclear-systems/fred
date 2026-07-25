#!/usr/bin/env python3
"""
Axial-max cladding wastage vs burnup for the two wastage models.

Reads the HDF5 outputs of the two sub-examples (run them first):
  precipitation_kinetics/fred_m_na_pk.h5
  la_tracking/fred_m_na_la.h5

Plots max-over-layers wastage [µm] against rod-average burnup [at%], with
the MFUEL operating limit delta = 0.5 * t_clad,0 as a reference line.

Run:  python plot_wastage_vs_burnup.py
"""

import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CASES = [
    ('PrecipitationKinetics (Clanth.for fit)',
     os.path.join(THIS_DIR, 'precipitation_kinetics', 'fred_m_na_pk.h5'), 'C0', '-'),
    ('LaTracking (MFUEL App. 9.3)',
     os.path.join(THIS_DIR, 'la_tracking', 'fred_m_na_la.h5'), 'C1', '--'),
]

D2S = 86400.0
fig, (ax_t, ax_b) = plt.subplots(1, 2, figsize=(12.5, 5))
for lbl, h5, color, ls in CASES:
    if not os.path.exists(h5):
        print(f"skipping {lbl}: {h5} not found (run the sub-example first)")
        continue
    with h5py.File(h5, 'r') as f:
        t  = f["time"][:]
        bu = f["burnup/bup_atpct"][:]
        xw = f["burnup/xwast_layer"][:]
    xmax = xw.max(axis=1) * 1e6
    ax_t.plot(t / D2S, xmax, ls, color=color, lw=1.6, label=lbl)
    ax_b.plot(bu, xmax, ls, color=color, lw=1.6, label=lbl)
    print(f"{lbl}: axial-max wastage at EOL = {xmax[-1]:.1f} µm")

ax_t.set_xlabel('time [days]')
ax_t.set_ylabel('axial-max cladding wastage [µm]')
ax_t.set_title('vs time (cf. Timpano Fig. 8.7: ~43 µm at 2176 d)')
ax_t.grid(alpha=0.3); ax_t.legend(loc='upper left')

ax_b.set_xlabel('rod-average burnup [at%]')
ax_b.set_ylabel('axial-max cladding wastage [µm]')
ax_b.set_title('vs burnup')
ax_b.grid(alpha=0.3); ax_b.legend(loc='upper left')

fig.suptitle('Cladding wastage: precipitation kinetics vs lanthanide tracking '
             '(Timpano benchmark specifications)')
fig.tight_layout()
out = os.path.join(THIS_DIR, 'wastage_vs_burnup.png')
fig.savefig(out, dpi=140)
print(f"Plot written: {out}")
