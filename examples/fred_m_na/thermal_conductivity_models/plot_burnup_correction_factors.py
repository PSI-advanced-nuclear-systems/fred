#!/usr/bin/env python3
"""
Burnup correction factors for U-Pu-Zr thermal conductivity.

Plots the multiplicative factor applied to k_fresh as a function of burnup
B [at%] for the two burnup-only models implemented in
src/apps/fred_m_na/fuelpelletmaterial/UPuZr.cpp:

  - EmpiricalBurnup : piecewise Karahan (2009) correction   (Flamb.for f=2)
  - SigmoidBurnup   : ESFR-SIMPLE double-logistic fit        (Flamb.for f=3)

Both factors are ~1 at B = 0 and stay within (0, 1].  The detailed
Na-infiltration model (f=1) is not shown: its correction depends on the
porosity/sodium state, not on burnup alone.

Run:  python plot_burnup_correction_factors.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def factor_piecewise(B_atpct):
    """Karahan (2009) piecewise burnup correction (UPuZr.cpp, EmpiricalBurnup)."""
    B = np.asarray(B_atpct, dtype=float)
    p = 0.135 * B
    f = np.where(B <= 2.0, (1.0 - p) / (1.0 + 1.7 * p),
        np.where(B <= 5.0, 0.5 + 0.0667 * (B - 2.0),
                 0.7))
    return f


def factor_sigmoid(B_atpct):
    """ESFR-SIMPLE double-logistic fit (UPuZr.cpp, SigmoidBurnup)."""
    B = np.asarray(B_atpct, dtype=float)
    a1, b1, c1 = 100.27110592, 0.5474772, -8.15852445
    a2, b2, c2 = 99.62701, 0.64761187, -6.4644582
    return (a1 / (1.0 + np.exp(-b1 * (B - c1)))
            - a2 / (1.0 + np.exp(-b2 * (B - c2))))


B = np.linspace(0.0, 20.0, 801)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(B, factor_piecewise(B), '-',
        label='EmpiricalBurnup — piecewise (Karahan 2009)')
ax.plot(B, factor_sigmoid(B), '--',
        label='SigmoidBurnup — double logistic (ESFR-SIMPLE)')

# Saturation levels and piecewise breakpoints
ax.axhline(0.7,   color='C0', ls=':', lw=0.8, alpha=0.6)
ax.axhline(0.644, color='C1', ls=':', lw=0.8, alpha=0.6)
ax.text(19.8, 0.706, '0.70', color='C0', ha='right', fontsize=9)
ax.text(19.8, 0.650, '0.644', color='C1', ha='right', fontsize=9)
for b in (2.0, 5.0):
    ax.axvline(b, color='C0', ls=':', lw=0.8, alpha=0.6)

ax.set_xlabel('burnup $B$ [at%]')
ax.set_ylabel('correction factor  $k_{corr}/k_{fresh}$  [-]')
ax.set_xlim(0.0, 20.0)
ax.grid(alpha=0.3)
ax.legend()
ax.set_title('U-Pu-Zr burnup correction to fresh-fuel thermal conductivity')

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'burnup_correction_factors.png')
fig.savefig(out, dpi=140)
print(f"Plot written: {out}")

print(f"\n{'B [at%]':>8} {'piecewise':>10} {'sigmoid':>10}")
for b in [0, 1, 2, 3, 5, 8, 12, 16, 20]:
    print(f"{b:8.1f} {factor_piecewise(b):10.3f} {factor_sigmoid(b):10.3f}")
