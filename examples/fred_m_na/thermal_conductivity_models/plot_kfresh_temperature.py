#!/usr/bin/env python3
"""
Fresh U-Pu-Zr thermal conductivity vs temperature (Aydin/Karahan 2009).

k_fresh(T) = A + B*T + C*T^2   [W/(m·K)], T in K, with composition-dependent
coefficients (upuzrFreshConductivity in UPuZr.cpp / Flamb.for):

  A = 17.5 * ((1 - 2.23*x_Zr)/(1 + 1.61*x_Zr) - 2.62*x_Pu)
  B = 1.54e-2 * (1 + 0.061*x_Zr)/(1 + 1.61*x_Zr)
  C = 9.38e-6 * (1 - 2.7*x_Pu)

Composition from the ESFR-SIMPLE example (x_Pu = 0.1325, x_Zr = 0.10 weight
fractions), plus +/-20% sensitivities on x_Pu (at fixed x_Zr) and on x_Zr
(at fixed x_Pu).  The x_Zr sensitivity matters for interpreting irradiated
radial k profiles: Zr redistribution enriches the hot gamma-phase centerline
and depletes the mid-radius, so the local k_fresh shifts along these curves.

Run:  python plot_kfresh_temperature.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def k_fresh(T_K, x_pu, x_zr):
    """Aydin/Karahan (2009) fresh U-Pu-Zr conductivity [W/(m·K)]."""
    A = 17.5 * ((1.0 - 2.23 * x_zr) / (1.0 + 1.61 * x_zr) - 2.62 * x_pu)
    B = 1.54e-2 * (1.0 + 0.061 * x_zr) / (1.0 + 1.61 * x_zr)
    C = 9.38e-6 * (1.0 - 2.7 * x_pu)
    return A + B * T_K + C * T_K**2


X_ZR = 0.10      # ESFR-SIMPLE example
X_PU = 0.1325

T = np.linspace(400.0, 1400.0, 501)   # up to ~solidus for this composition

fig, ax = plt.subplots(figsize=(8.5, 5.5))
ax.plot(T, k_fresh(T, X_PU, X_ZR), 'k-', lw=2,
        label=f'nominal  ($x_{{Pu}}$={X_PU:.4f}, $x_{{Zr}}$={X_ZR:.2f})')
ax.plot(T, k_fresh(T, 0.8 * X_PU, X_ZR), '--',  color='C0',
        label=f'$x_{{Pu}}$ -20%  ({0.8*X_PU:.4f})')
ax.plot(T, k_fresh(T, 1.2 * X_PU, X_ZR), '-.',  color='C0',
        label=f'$x_{{Pu}}$ +20%  ({1.2*X_PU:.4f})')
ax.plot(T, k_fresh(T, X_PU, 0.8 * X_ZR), '--',  color='C1',
        label=f'$x_{{Zr}}$ -20%  ({0.8*X_ZR:.3f})')
ax.plot(T, k_fresh(T, X_PU, 1.2 * X_ZR), '-.',  color='C1',
        label=f'$x_{{Zr}}$ +20%  ({1.2*X_ZR:.3f})')

ax.set_xlabel('temperature [K]')
ax.set_ylabel('$k_{fresh}$  [W/(m·K)]')
ax.set_title('Fresh U-Pu-Zr conductivity (Aydin/Karahan 2009), '
             'composition sensitivity')
ax.grid(alpha=0.3)
ax.legend()

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'kfresh_vs_temperature.png')
fig.savefig(out, dpi=140)
print(f"Plot written: {out}")

print(f"\n{'T [K]':>7} {'nominal':>9} {'-20%Pu':>9} {'+20%Pu':>9} "
      f"{'-20%Zr':>9} {'+20%Zr':>9}")
for t in [500, 700, 900, 1100, 1300]:
    print(f"{t:7.0f} {k_fresh(t, X_PU, X_ZR):9.2f} "
          f"{k_fresh(t, 0.8*X_PU, X_ZR):9.2f} {k_fresh(t, 1.2*X_PU, X_ZR):9.2f} "
          f"{k_fresh(t, X_PU, 0.8*X_ZR):9.2f} {k_fresh(t, X_PU, 1.2*X_ZR):9.2f}")
