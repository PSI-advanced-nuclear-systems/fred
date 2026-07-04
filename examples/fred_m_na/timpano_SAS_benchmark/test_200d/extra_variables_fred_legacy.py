"""
Parse legacy FRED-M output files (rstfrd*) from the IF-AVG 1000s benchmark
and extract additional variables for comparison with FRED-M-Na:

  - xwast      : lanthanide cladding wastage (m and µm), per axial node
  - hoop_strain: total cladding hoop strain (%), per axial × radial node
  - fgr_moles  : fission gas release (mol), per axial node × radial FGR zone
                 The 2-D shape (nz=24, n_fgr_zones=10) arises because the
                 FEAST/GRSIS model tracks release separately in each radial
                 fuel ring; summing over axis=1 gives total moles per slice.
  - gpres      : plenum gas pressure (Pa and MPa), single scalar

Usage:
    from extra_variables_fred_legacy import snapshots
    # snapshots is an ordered dict: filename -> dict of arrays/scalars

    # Example: peak wastage vs time
    import numpy as np
    times = [snapshots[k]["time_d"] for k in sorted(snapshots)]
    xwast_peak = [snapshots[k]["xwast_um"].max() for k in sorted(snapshots)]
"""

import os
import re
import numpy as np

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "legacy")

_NZ = 24   # number of axial nodes
_NR_HOOP = 16   # fuel(11) + cladding(5) radial nodes for hoop strain
_NR_FGR  = 10   # radial FGR zones used in FEAST/GRSIS model


def _parse_rstfrd(path: str) -> dict:
    with open(path) as fh:
        lines = fh.readlines()

    result: dict = {}

    # Time
    m = re.match(r"time\s+([\d.E+\-]+)\s+s\s+\|\s+([\d.E+\-]+)\s+d", lines[0])
    result["time_s"] = float(m.group(1))
    result["time_d"] = float(m.group(2))

    # Locate section header lines
    hoop_start = fgr_start = summ_start = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped == "total hoop strain (%)":
            hoop_start = i + 1
        elif stripped.startswith("fission gas release (moles)"):
            fgr_start = i + 1
        elif stripped.startswith("# ") and "xwast" in stripped:
            summ_start = i + 1

    # --- hoop strain: nz x nf+nc ---
    hoop = []
    for i in range(hoop_start, hoop_start + _NZ):
        # first token is the row index; rest are values
        vals = [float(x) for x in lines[i].split()[1:]]
        hoop.append(vals)
    result["hoop_strain_pct"] = np.array(hoop)   # shape (24, 16)

    # --- fission gas release (moles): nz x n_fgr_zones ---
    # Each row is one axial slice; columns are the 10 radial FGR zones.
    # The 2-D layout is because the FEAST model tracks gas per radial ring;
    # summing over zones gives total moles released per axial slice.
    fgr = []
    for i in range(fgr_start, fgr_start + _NZ):
        vals = [float(x) for x in lines[i].split()[1:]]
        fgr.append(vals)
    result["fgr_moles"] = np.array(fgr)   # shape (24, 10)

    # --- summary table: xwast, fgr_frac, gpres ---
    # Column mapping (0-based after stripping the leading row number):
    # 0=zzz  1=tfin  2=tfout  3=tcin  4=tcout  5=tfave  6=tcave  7=tcool
    # 8=ma   9=ntot  10=clfue 11=xwast 12=dose 13=neufl
    # 14=qs  15=qv   16=ql   17=ql2   18=htc  19=hgap  20=bup   21=bup2
    # 22=fggen 23=fgr  24=gap  25=gapth 26=pfc 27=vel  28=rfi
    # 29=rfo  30=rci  31=rco  32=rof0  33=ecsw 34=gask
    # 35=hgap1 36=hgap2 37=hgap3 38=ajump 39=dzf 40=dzc
    # 41=gpres 42=reloc 43=rloc 44=state
    #
    # parts[0] = row index, so parts[col+1] = column value
    #
    # NOTE: cols 22/23 ("fggen"/"fgr") are the rod-scalar cumulative gas
    # generated/released in MOLES (fggen(l)/fgrel(l) in Rstfrd.for), not a
    # percentage — the same value is repeated on every axial row. The true
    # FGR fraction is fgrel(l)/fggen(l); it must NOT be used directly as a
    # "percent" (multiplying the raw mole count by 100 previously produced
    # a ~13x-too-low number that was mistaken for a physics bug).
    xwast_vals  = []
    fggen_vals  = []
    fgrel_vals  = []
    gpres_vals  = []
    for i in range(summ_start, summ_start + _NZ):
        parts = lines[i].split()
        xwast_vals.append(float(parts[12]))    # col 11 → parts[12]
        fggen_vals.append(float(parts[23]))    # col 22 → parts[23], moles
        fgrel_vals.append(float(parts[24]))    # col 23 → parts[24], moles
        gpres_vals.append(float(parts[42]))    # col 41 → parts[42]

    result["xwast_m"]    = np.array(xwast_vals)          # (24,) m
    result["xwast_um"]   = np.array(xwast_vals) * 1e6    # (24,) µm
    result["fggen_mol"]  = np.array(fggen_vals)          # (24,) mol (rod-scalar, repeated)
    result["fgrel_mol"]  = np.array(fgrel_vals)          # (24,) mol (rod-scalar, repeated)
    result["fgr_frac"]   = result["fgrel_mol"] / result["fggen_mol"]  # (24,) true fraction
    # gpres is uniform across all axial rows (single plenum)
    result["gpres_Pa"]  = float(gpres_vals[0])
    result["gpres_MPa"] = float(gpres_vals[0]) * 1e-6

    return result


def _load_all() -> dict:
    files = sorted(
        f for f in os.listdir(_DATA_DIR) if f.startswith("rstfrd")
    )
    out = {}
    for fname in files:
        out[fname] = _parse_rstfrd(os.path.join(_DATA_DIR, fname))
    return out


# Module-level dict loaded at import time.
# Keys  : rstfrd file names (e.g. "rstfrd000000000021")
# Values: dict with keys time_s, time_d, xwast_m, xwast_um,
#         hoop_strain_pct, fgr_moles, fgr_frac, gpres_Pa, gpres_MPa
snapshots: dict = _load_all()


if __name__ == "__main__":
    # Quick sanity-check printout
    header = (
        f"{'file':<30} {'t (d)':>10} {'xwast_peak (µm)':>18} "
        f"{'gpres (MPa)':>14} {'FGR frac max (%)':>18}"
    )
    print(header)
    print("-" * len(header))
    for key in sorted(snapshots):
        d = snapshots[key]
        print(
            f"{key:<30} {d['time_d']:>10.1f} {d['xwast_um'].max():>18.1f} "
            f"{d['gpres_MPa']:>14.3f} {d['fgr_frac'].max()*100:>18.2f}"
        )
