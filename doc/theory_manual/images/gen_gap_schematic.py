"""Generate gap model schematic for the FRED theory manual."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fig, ax = plt.subplots(figsize=(8, 4))

# Fuel outer surface (left) and clad inner surface (right)
r_fo = 1.0
r_ci = 2.2
gap  = r_ci - r_fo

z_lo, z_hi = 0.0, 3.0

# Fill gap region with light blue (gas)
rect_gas = plt.Rectangle((r_fo, z_lo), gap, z_hi - z_lo,
                          facecolor='#e8f4fd', edgecolor='none', zorder=1)
ax.add_patch(rect_gas)

# Fuel region (left)
rect_fuel = plt.Rectangle((0.0, z_lo), r_fo, z_hi - z_lo,
                            facecolor='#d4e8ff', edgecolor='#336699', linewidth=1.2, zorder=2)
ax.add_patch(rect_fuel)
ax.text(0.5*r_fo, 0.5*(z_lo+z_hi), 'Fuel\n(outer)', ha='center', va='center', fontsize=9)

# Clad region (right)
rect_clad = plt.Rectangle((r_ci, z_lo), 1.5, z_hi - z_lo,
                            facecolor='#ffe0b2', edgecolor='#8B4513', linewidth=1.2, zorder=2)
ax.add_patch(rect_clad)
ax.text(r_ci + 0.75, 0.5*(z_lo+z_hi), 'Cladding\n(inner)', ha='center', va='center', fontsize=9)

# Roughness bumps on fuel outer surface
np.random.seed(42)
z_rough = np.linspace(z_lo + 0.2, z_hi - 0.2, 20)
fuel_rough_amp = 0.04
clad_rough_amp = 0.03
for zz in z_rough:
    bump_f = fuel_rough_amp * (0.5 + np.random.rand())
    bump_c = clad_rough_amp * (0.5 + np.random.rand())
    ax.plot(r_fo + bump_f, zz, '.', color='#336699', ms=2, zorder=4)
    ax.plot(r_ci - bump_c, zz, '.', color='#8B4513', ms=2, zorder=4)

# Gap width arrow
zmid = 0.5 * (z_lo + z_hi)
ax.annotate('', xy=(r_ci, zmid), xytext=(r_fo, zmid),
            arrowprops=dict(arrowstyle='<->', color='#222', lw=1.2))
ax.text(0.5*(r_fo + r_ci), zmid + 0.12, r'$\delta_\mathrm{gap}$',
        ha='center', va='bottom', fontsize=10)

# Roughness annotation
z_ruf = z_hi - 0.5
ax.annotate('', xy=(r_fo + 0.06, z_ruf), xytext=(r_fo, z_ruf),
            arrowprops=dict(arrowstyle='<->', color='#336699', lw=0.8))
ax.text(r_fo + 0.08, z_ruf, r'$\rho_f$', ha='left', va='center', fontsize=9, color='#336699')

ax.annotate('', xy=(r_ci - 0.04, z_ruf + 0.2), xytext=(r_ci, z_ruf + 0.2),
            arrowprops=dict(arrowstyle='<->', color='#8B4513', lw=0.8))
ax.text(r_ci - 0.06, z_ruf + 0.2, r'$\rho_c$', ha='right', va='center', fontsize=9, color='#8B4513')

# Temperature labels
ax.text(r_fo - 0.08, z_hi - 0.2, r'$T_f$', ha='right', va='center',
        fontsize=10, color='#003366', fontweight='bold')
ax.text(r_ci + 0.08, z_hi - 0.2, r'$T_c$', ha='left', va='center',
        fontsize=10, color='#6B2600', fontweight='bold')

# Gas label
ax.text(0.5*(r_fo+r_ci), 0.18*(z_lo+z_hi), 'Fill gas\n(He + FG)',
        ha='center', va='center', fontsize=8.5, color='#1565c0',
        fontstyle='italic')

# Heat flux arrow
ax.annotate('', xy=(r_ci + 0.2, 1.5), xytext=(r_fo - 0.2, 1.5),
            arrowprops=dict(arrowstyle='->', color='red', lw=1.5))
ax.text(0.5*(r_fo + r_ci), 1.5 - 0.15, r"$q'' = h_\mathrm{gap}(T_f - T_c)$",
        ha='center', va='top', fontsize=8.5, color='red')

ax.set_xlim(-0.15, r_ci + 1.65)
ax.set_ylim(z_lo - 0.1, z_hi + 0.25)
ax.set_xlabel(r'Radial position $r$', fontsize=10)
ax.set_ylabel(r'Axial position $z$', fontsize=10)
ax.set_title('Fuel--Cladding Gap Model', fontsize=10)
ax.grid(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('gap_schematic.pdf', bbox_inches='tight', dpi=150)
print("Saved gap_schematic.pdf")
