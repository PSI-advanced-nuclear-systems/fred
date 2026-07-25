#!/usr/bin/env python3
"""
Compare cladding void swelling: Hofmann vs SAS4A (platform), overlaid with
the legacy FRED-M reference (Cswel.for icswel=1/icswel=2) at the same
IF-AVG axial layer.

Legacy reference paths default to the scratchpad verification decks used
to validate this port; override LEGACY_DIR if you re-ran them elsewhere.
Legacy runs were 400 d (shorter, for fast turnaround); platform runs are
the full 2176 d benchmark — the comparison plot is clipped to the overlap.

Run:  python compare_clad_swelling.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from common import read_swelling, load_if_avg_deck, D2S

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
LEGACY_DIR = os.environ.get(
    'LEGACY_DIR',
    '/tmp/claude-1000/-mnt-c-Users-chan-y-Documents-fred-platform-'
    'examples-fred-m-na-thermal-conductivity-models/'
    '2857ab32-3c6d-44d3-9d63-6f15897248e8/scratchpad/clad_swelling_verify')

deck = load_if_avg_deck()
NZ, NC = deck['nz'], deck['nc']
LAYER = 11   # legacy OUTPUT_AXIAL_LAYER 12 -> 0-based 11
NODE  = 0    # legacy ecs(1,j,l) -> innermost clad interval, 0-based 0


def read_legacy(subdir):
    fd = os.path.join(LEGACY_DIR, subdir, 'fred.dat')
    if not os.path.exists(fd):
        return None
    with open(fd) as f:
        hdr = f.readline().split()
    cols = {n: i for i, n in enumerate(hdr)}
    raw = np.genfromtxt(fd, skip_header=1, dtype=str)
    t   = raw[:, cols['time(s)']].astype(float)
    ecs = raw[:, cols['ecs(%)']].astype(float) / 100.0   # -> strain [-]
    dose = raw[:, cols['dpa(dpa)']].astype(float)
    return t, ecs, dose


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

for model_name, color in [('hofmann', 'C0'), ('sas4a', 'C1')]:
    h5 = os.path.join(THIS_DIR, model_name, f'{model_name}.h5')
    if not os.path.exists(h5):
        print(f"skipping {model_name}: {h5} not found (run the sub-example first)")
        continue
    t, ecs, dose, neuflue2 = read_swelling(h5, NZ, NC, LAYER, NODE)
    ax1.plot(t / D2S, ecs * 100.0, '-', color=color, marker='o',
             markevery=max(len(t)//20, 1), markersize=5, markerfacecolor='none',
             label=f'platform {model_name}')
    print(f"platform {model_name}: ecs @ EOL = {ecs[-1]*100:.4f} %  "
          f"dose @ EOL = {dose[-1]:.1f} dpa")

    leg = read_legacy(model_name)
    if leg is not None:
        lt, lecs, ldose = leg
        ax1.plot(lt / D2S, lecs * 100.0, '--', color=color, marker='x',
                 markevery=max(len(lt)//20, 1), markersize=6,
                 label=f'legacy {model_name}')
        print(f"legacy   {model_name}: ecs @ {lt[-1]/D2S:.0f} d = "
              f"{lecs[-1]*100:.4f} %  dose = {ldose[-1]:.1f} dpa")
        # overlay dose on right panel too
        ax2.plot(lt / D2S, ldose, '--', color=color, marker='x',
                 markevery=max(len(lt)//20, 1), markersize=6,
                 label=f'legacy {model_name}')

    ax2.plot(t / D2S, dose, '-', color=color, marker='o',
             markevery=max(len(t)//20, 1), markersize=5, markerfacecolor='none',
             label=f'platform {model_name}')

ax2.axhline(100.0, color='gray', ls=':', lw=1.0)
ax2.text(5, 105, 'SAS4A onset (100 dpa)', fontsize=8, color='gray')

ax1.set_xlabel('time [days]')
ax1.set_ylabel('cladding void swelling strain, ecs [%]')
ax1.set_title(f'Void swelling (IF-AVG, layer {LAYER+1}, clad node {NODE+1})')
ax1.grid(alpha=0.3); ax1.legend(fontsize=8)

ax2.set_xlabel('time [days]')
ax2.set_ylabel('cladding dose [dpa]')
ax2.set_title('Cladding dose (shared SAS fluence/dose bookkeeping)')
ax2.grid(alpha=0.3); ax2.legend(fontsize=8)

fig.suptitle('Cladding void swelling: Hofmann vs SAS4A, platform vs legacy')
fig.tight_layout()
out = os.path.join(THIS_DIR, 'clad_swelling_comparison.png')
fig.savefig(out, dpi=140)
print(f"\nPlot written: {out}")
