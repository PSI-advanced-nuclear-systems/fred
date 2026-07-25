#!/usr/bin/env python3
"""
Mesh-refinement study for the LaTracking cladding-wastage model
(theory manual figure images/la_mesh_refinement.pdf).

Standalone replica of cladWastageLaTracking (FredMNaCladWastage.hpp) driven
by the IF-AVG layer-12 temperature profile at the Timpano benchmark power:
wastage at 1825 d vs radial subgrid refinement of the nf=11 thermal grid,
plus the end-of-life near-interface concentration profiles showing the
boundary layer that the coarse grid cannot resolve.

Run from doc/theory_manual/images/:  python make_la_mesh_refinement.py
(requires the lanthanide_wastage benchmark H5 for the temperature profile;
falls back to a representative built-in profile if absent)
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RFO = 3.863e-3
NF  = 11
Y_LA, RG   = 15.29e9, 8.314
D0, Q, TFL = 5.0, 250.0e3, 800.0
U_SOL      = 6.4e27
QV         = 1.022e9
TEND, DT   = 1825 * 86400.0, 86400.0

# Layer-12 fuel temperature profile (IF-AVG, mid-life) — from the benchmark
# H5 if available, otherwise the recorded values.
T11 = np.array([1009.6, 1006.6, 997.8, 983.4, 963.7, 939.2,
                910.5, 878.2, 843.1, 805.6, 784.2])
h5 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  '..', '..', '..', 'examples', 'fred_m_na',
                  'timpano_SAS_benchmark', 'lanthanide_wastage',
                  'if_avg_la.h5')
try:
    import h5py
    with h5py.File(h5, 'r') as f:
        t = f['time'][:]
        T = f['thermal/T'][:]
    T11 = T[np.argmin(abs(t - 900 * 86400)), 11 * 16:11 * 16 + NF]
    print("temperature profile taken from", os.path.basename(h5))
except Exception as e:
    print("using built-in T profile:", e)


def run(refine):
    n  = (NF - 1) * refine + 1
    dr = RFO / (n - 1)
    r  = np.linspace(0.0, RFO, n)
    Tn = np.interp(r, np.linspace(0, RFO, NF), T11)
    D  = D0 * np.exp(-Q / (RG * np.maximum(Tn, TFL)))
    V  = np.zeros(n)
    V[0] = 0.5 * (0.5 * dr) ** 2
    for k in range(1, n - 1):
        V[k] = (k * dr) * dr
    V[-1] = ((n - 1) * dr - 0.25 * dr) * (0.5 * dr)
    S  = Y_LA * QV
    h_stab = 0.25 * dr * dr / (4 * D.max())
    rf = (np.arange(n - 1) + 0.5) * dr
    Df = 0.5 * (D[:-1] + D[1:])
    c  = np.zeros(n)
    xw = 0.0
    t_ = 0.0
    while t_ < TEND:
        h = min(DT, h_stab, TEND - t_)
        flux = rf * Df * (c[:-1] - c[1:]) / dr
        sink = (n - 1) * dr * D[-1] * c[-1] / (0.5 * dr)
        inv0 = (V * c).sum()
        cn = c.copy()
        cn[0]    += h * (-flux[0] / V[0] + S)
        cn[1:-1] += h * ((flux[:-1] - flux[1:]) / V[1:-1] + S)
        cn[-1]   += h * ((flux[-1] - sink) / V[-1] + S)
        c = np.maximum(cn, 0)
        absorbed = S * h * V.sum() - ((V * c).sum() - inv0)
        if absorbed > 0:
            xw += absorbed / (RFO * U_SOL)
        t_ += h
    return xw * 1e6, r, c


refines = [1, 2, 4, 8, 16]
results = {m: run(m) for m in refines}
xw16 = results[16][0]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

xws = [results[m][0] for m in refines]
ax1.plot(refines, xws, 'o-')
for m, w in zip(refines, xws):
    ax1.annotate(f'{w:.1f}', (m, w), textcoords='offset points',
                 xytext=(6, -10), fontsize=8)
ax1.axvline(8, color='C1', ls='--', lw=1.0)
ax1.text(8 * 1.05, xws[0] + 2, 'default\nla_refine = 8', fontsize=8, color='C1')
ax1.set_xscale('log', base=2)
ax1.set_xticks(refines, [str(m) for m in refines])
ax1.set_xlabel('subgrid refinement factor per thermal cell')
ax1.set_ylabel('wastage at 1825 d [µm]')
ax1.set_title('(a) Mesh convergence (IF-AVG, layer 12)')
ax1.grid(alpha=0.3)

for m, style in [(1, 'o-'), (8, 's-'), (16, '.-')]:
    _, r, c = results[m]
    ax2.plot(r * 1e3, c / 1e27, style, ms=3,
             label=f'refine {m} ({(NF-1)*m+1} pts)')
ax2.set_xlim(3.0, RFO * 1e3)
ax2.set_xlabel('radius [mm]')
ax2.set_ylabel('La density $u$ [$10^{27}$ #/m³]')
ax2.set_title('(b) EOL near-interface boundary layer')
ax2.grid(alpha=0.3)
ax2.legend(fontsize=8)

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'la_mesh_refinement.pdf')
fig.savefig(out)
print("wrote", out)
print({m: round(results[m][0], 1) for m in refines},
      " rel-to-16:", [f"{results[m][0]/xw16:.3f}" for m in refines])
