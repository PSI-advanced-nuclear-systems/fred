"""Generate r-z discretization schematic for the FRED theory manual."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(figsize=(9, 6))

# Geometry parameters (scaled for illustration, nf=5, nc=3, nz=3)
nf = 5
nc = 3
nz = 3

r_fi = 0.0
r_fo = 4.2
r_ci = 4.5  # gap
r_co = 5.4
dz   = 2.0

# Draw axial layers
z_edges = [j * dz for j in range(nz + 1)]

colors_fuel = ['#d4e8ff', '#c0dbff', '#abd0ff', '#92c2ff', '#78b3ff']
colors_clad = ['#ffe0b2', '#ffcc88', '#ffb860']

for j in range(nz):
    zlo, zhi = z_edges[j], z_edges[j+1]
    # Fuel rings
    dr_f = (r_fo - r_fi) / nf
    for i in range(nf):
        r_lo = r_fi + i * dr_f
        r_hi = r_fi + (i+1) * dr_f
        rect = plt.Rectangle((r_lo, zlo), r_hi - r_lo, dz,
                              facecolor=colors_fuel[i], edgecolor='#336699', linewidth=0.6)
        ax.add_patch(rect)
        # Node dot at centre
        ax.plot(0.5*(r_lo+r_hi), 0.5*(zlo+zhi), 'o', color='#003366', ms=4, zorder=5)

    # Gap region
    rect_gap = plt.Rectangle((r_fo, zlo), r_ci - r_fo, dz,
                               facecolor='#f0f0f0', edgecolor='#888', linewidth=0.6,
                               linestyle='--')
    ax.add_patch(rect_gap)
    ax.text(0.5*(r_fo+r_ci), 0.5*(zlo+zhi), 'gap', ha='center', va='center',
            fontsize=7, color='#555')

    # Cladding rings
    dr_c = (r_co - r_ci) / nc
    for i in range(nc):
        r_lo = r_ci + i * dr_c
        r_hi = r_ci + (i+1) * dr_c
        rect = plt.Rectangle((r_lo, zlo), r_hi - r_lo, dz,
                              facecolor=colors_clad[i], edgecolor='#8B4513', linewidth=0.6)
        ax.add_patch(rect)
        ax.plot(0.5*(r_lo+r_hi), 0.5*(zlo+zhi), 's', color='#6B2600', ms=4, zorder=5)

# Radius annotations (top row only)
zref = z_edges[nz] + 0.2
fontsize_ann = 8.5

def vline(r, label, yb, yt, color='black'):
    ax.annotate('', xy=(r, yt), xytext=(r, yb),
                arrowprops=dict(arrowstyle='->', color=color, lw=0.8))
    ax.text(r, yt + 0.08, label, ha='center', va='bottom', fontsize=fontsize_ann, color=color)

ax.axhline(z_edges[nz], color='#aaa', linestyle=':', linewidth=0.5)
for r, lbl in [(0.0, r'$r_{fi}$'), (r_fo, r'$r_{fo}$'),
               (r_ci, r'$r_{ci}$'), (r_co, r'$r_{co}$')]:
    ax.plot([r, r], [0, z_edges[nz]], color='gray', linestyle=':', linewidth=0.8, zorder=2)
    ax.text(r, z_edges[nz] + 0.25, lbl, ha='center', va='bottom', fontsize=9)

# Axial layer labels
for j in range(nz):
    zmid = 0.5*(z_edges[j] + z_edges[j+1])
    ax.text(-0.5, zmid, rf'$j={j}$', ha='right', va='center', fontsize=8, color='#222')
    ax.annotate('', xy=(r_co + 0.25, z_edges[j+1]), xytext=(r_co + 0.25, z_edges[j]),
                arrowprops=dict(arrowstyle='<->', color='#444', lw=0.8))
    ax.text(r_co + 0.55, zmid, r'$\Delta z_j$', ha='left', va='center', fontsize=8)

# Fuel ring labels (first axial layer)
dr_f = (r_fo - r_fi) / nf
j = 0
zlo = z_edges[j]
for i in range(nf):
    r_mid = r_fi + (i + 0.5) * dr_f
    ax.text(r_mid, zlo - 0.18, rf'$i={i}$', ha='center', va='top', fontsize=7, color='#003366')

dr_c = (r_co - r_ci) / nc
for i in range(nc):
    r_mid = r_ci + (i + 0.5) * dr_c
    ax.text(r_mid, zlo - 0.18, rf'$i={i}$', ha='center', va='top', fontsize=7, color='#6B2600')

# Legend
fuel_patch = mpatches.Patch(facecolor=colors_fuel[2], edgecolor='#336699', label=f'Fuel rings ($n_f={nf}$)')
clad_patch = mpatches.Patch(facecolor=colors_clad[1], edgecolor='#8B4513', label=f'Cladding rings ($n_c={nc}$)')
gap_patch  = mpatches.Patch(facecolor='#f0f0f0', edgecolor='#888', linestyle='--', label='Gap')
ax.legend(handles=[fuel_patch, clad_patch, gap_patch],
          loc='upper right', fontsize=8, framealpha=0.9)

ax.set_xlabel(r'Radial position $r$ [mm]', fontsize=10)
ax.set_ylabel(r'Axial position $z$ [arb.]', fontsize=10)
ax.set_title(r'FRED Discretization: $n_f$ fuel rings, $n_c$ cladding rings, $n_z$ axial layers',
             fontsize=10)
ax.set_xlim(-1.0, r_co + 1.4)
ax.set_ylim(-0.6, z_edges[nz] + 0.8)
ax.set_aspect('equal')
ax.grid(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('rz_discretization.pdf', bbox_inches='tight', dpi=150)
print("Saved rz_discretization.pdf")
