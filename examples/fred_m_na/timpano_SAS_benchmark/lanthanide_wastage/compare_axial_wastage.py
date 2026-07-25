#!/usr/bin/env python3
"""
Axial lanthanide-wastage profiles at 1825 days for the three Timpano
benchmark cases (cf. Timpano's FRED-MFUEL comparison figure, top row):

  OF-AVG : dome ~12-28 um            IF-AVG : dome ~30-46 um
  IF-PEAK: rising to ~100-110 um at the top of the core

Reads the H5 pairs written by run_if_avg.py / run_if_peak.py / run_of_avg.py
and plots PrecipitationKinetics vs LaTracking per case.

Run:  python compare_axial_wastage.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from common import axial_wastage_at, D2S

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
T_CMP    = 1825 * D2S

CASES = [('of_avg', 'OF-AVG'), ('if_avg', 'IF-AVG'), ('if_peak', 'IF-PEAK')]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
for ax, (case, title) in zip(axes, CASES):
    for model, lbl, style in [('pk', 'PrecipitationKinetics', 'o-'),
                              ('la', 'LaTracking (MFUEL)',    's--')]:
        h5 = os.path.join(THIS_DIR, f'{case}_{model}.h5')
        if not os.path.exists(h5):
            print(f"skipping {case}/{model}: {h5} not found")
            continue
        try:
            nodes, xw, t_actual = axial_wastage_at(h5, T_CMP)
        except OSError as e:
            print(f"skipping {case}/{model}: file locked/incomplete ({e})")
            continue
        if t_actual < 0.95 * T_CMP:
            print(f"skipping {case}/{model}: run only reached {t_actual/D2S:.0f} d")
            continue
        ax.plot(nodes, xw * 1e6, style, ms=4, label=lbl)
        print(f"{title:8s} {lbl:22s} @ {t_actual/D2S:6.0f} d: "
              f"min {xw.min()*1e6:6.1f}  max {xw.max()*1e6:6.1f} µm")
    ax.set_xlabel('axial node')
    ax.set_ylabel('lanthanide wastage [µm]')
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

fig.suptitle('Axial lanthanide wastage at 1825 days '
             '(cf. Timpano FRED-MFUEL benchmark)')
fig.tight_layout()
out = os.path.join(THIS_DIR, 'axial_wastage_1825d.png')
fig.savefig(out, dpi=140)
print(f"Plot written: {out}")
