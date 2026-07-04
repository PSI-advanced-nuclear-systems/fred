#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    r = np.linspace(0.0, 1.0, 500)  # normalized radius r/R

    # Illustrative MOX radial temperature profile (degC): hot center, cooler rim.
    T_center = 2400.0
    T_surface = 750.0
    temp = T_surface + (T_center - T_surface) * (1.0 - r**1.8)

    # Illustrative zoning radii from a representative operating state.
    r_crsz = 0.30
    r_unres = 0.80
    T_restructure = 1800.0

    # Illustrative radial burnup profile (arb.): flatter center, rim increase.
    burnup = 0.65 + 0.35 * r**1.7

    # Illustrative zone-wise FGR components used only for conceptual plotting.
    fgr_crsz = 0.85 * np.exp(-((r / max(r_crsz, 1.0e-6)) ** 3.2))
    fgr_unres = 0.20 * (1.0 - np.exp(-3.0 * np.clip(r - r_crsz, 0.0, None)))
    fgr_unres *= 1.0 / (1.0 + np.exp(-(burnup - 0.82) / 0.04))
    fgr_total = np.clip(fgr_crsz + fgr_unres, 0.0, 1.0)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), constrained_layout=True)

    # Panel A: temperature and zones.
    ax = axes[0]
    ax.plot(r, temp, color="#1f4e79", lw=2.4, label="MOX temperature")
    ax.axhline(T_restructure, color="#7a7a7a", lw=1.5, ls="--", label="Restructuring threshold")

    ax.axvspan(0.0, r_crsz, color="#f4b183", alpha=0.35, label="CRSZ")
    ax.axvspan(r_crsz, r_unres, color="#c6e0b4", alpha=0.35, label="Unrestructured")
    ax.axvspan(r_unres, 1.0, color="#d9d2e9", alpha=0.35, label="As-fabricated")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(600.0, 2550.0)
    ax.set_xlabel("Normalized radius r/R")
    ax.set_ylabel("Temperature [degC]")
    ax.set_title("(a) Temperature profile and radial zones")

    ax.text(0.06, 2360.0, "CRSZ", fontsize=10, weight="bold")
    ax.text(0.43, 2360.0, "Unrestructured", fontsize=10, weight="bold")
    ax.text(0.84, 2360.0, "As-fabricated", fontsize=10, weight="bold")
    ax.legend(loc="lower left", fontsize=8, frameon=True)

    # Panel B: burnup and zone-wise FGR behavior.
    ax2 = axes[1]
    ax2.plot(r, burnup, color="#7f6000", lw=2.2, label="Typical radial burnup (arb.)")
    ax2.plot(r, fgr_crsz, color="#b22222", lw=2.0, label="FGR contribution: CRSZ")
    ax2.plot(r, fgr_unres, color="#2e8b57", lw=2.0, label="FGR contribution: unrestructured")
    ax2.plot(r, fgr_total, color="#111111", lw=2.5, ls="-.", label="Total FGR (illustrative)")

    ax2.axvline(r_crsz, color="#666666", ls=":", lw=1.3)
    ax2.axvline(r_unres, color="#666666", ls=":", lw=1.3)

    ax2.set_xlim(0.0, 1.0)
    ax2.set_ylim(0.0, 1.05)
    ax2.set_xlabel("Normalized radius r/R")
    ax2.set_ylabel("Normalized value [-]")
    ax2.set_title("(b) Typical burnup and zone-wise FGR behavior")
    ax2.legend(loc="upper right", fontsize=8, frameon=True)

    fig.suptitle("Waltar-Reynolds conceptual zoning for MOX fuel", fontsize=12, weight="bold")
    fig.savefig("waltar_reynolds_zones.pdf")
    print("Saved waltar_reynolds_zones.pdf")


if __name__ == "__main__":
    main()
