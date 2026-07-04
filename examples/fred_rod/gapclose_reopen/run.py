"""
Gap close / reopen — gas bond verification against legacy FRED.

Gas pressure evolves via the GasBondPressure DAE (ideal gas law with the
current gap volume and plenum temperature) so it rises with temperature and
drops back as the fuel contracts.

Scenario:
  Phase 1 (0-100 s):    q = 0         → open-gap steady state
  Phase 2 (100-200 s):  ramp 0→5e8    → gap closes (thermal expansion)
  Phase 2 (200-600 s):  q = 5e8 W/m3  → closed-gap steady state
  Phase 3 (600-700 s):  ramp 5e8→0    → power down
  Phase 3 (700-1500 s): q = 0         → fuel contracts, gap reopens

Geometry: solid dummy pellet, nf=3, nc=2, nz=1
  rfo0=4.2 mm, gap=15 µm (rci0=4.215 mm), rco0=4.9 mm, ruff=rufc=5 µm
"""

import sys, os, glob, re
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'build'))
import fred_rod as fred

# ---------------------------------------------------------------------------
# Embedded legacy reference data (fallback if outfrd files are unavailable)
# ---------------------------------------------------------------------------
LEGACY_TIMES = [
    0.0, 10.0, 20.0, 30.0, 40.0, 50.0,
    60.0, 70.0, 80.0, 90.0, 100.0, 110.0,
    120.0, 130.0, 140.0, 150.0, 160.0, 170.0,
    180.0, 190.0, 200.0, 210.0, 220.0, 230.0,
    240.0, 250.0, 260.0, 270.0, 280.0, 290.0,
    300.0, 310.0, 320.0, 330.0, 340.0, 350.0,
    360.0, 370.0, 380.0, 390.0, 400.0, 410.0,
    420.0, 430.0, 440.0, 450.0, 460.0, 470.0,
    480.0, 490.0, 500.0, 510.0, 520.0, 530.0,
    540.0, 550.0, 560.0, 570.0, 580.0, 590.0,
    600.0, 610.0, 620.0, 630.0, 640.0, 650.0,
    660.0, 670.0, 680.0, 690.0, 700.0, 710.0,
    720.0, 730.0, 740.0, 750.0, 760.0, 770.0,
    780.0, 790.0, 800.0, 810.0, 820.0, 830.0,
    840.0, 850.0, 860.0, 870.0, 880.0, 890.0,
    900.0, 910.0, 920.0, 930.0, 940.0, 950.0,
    960.0, 970.0, 980.0, 990.0, 1000.0, 1010.0,
    1020.0, 1030.0, 1040.0, 1050.0, 1060.0, 1070.0,
    1080.0, 1090.0, 1100.0, 1110.0, 1120.0, 1130.0,
    1140.0, 1150.0, 1160.0, 1170.0, 1180.0, 1190.0,
    1200.0, 1210.0, 1220.0, 1230.0, 1240.0, 1250.0,
    1260.0, 1270.0, 1280.0, 1290.0, 1300.0, 1310.0,
    1320.0, 1330.0, 1340.0, 1350.0, 1360.0, 1370.0,
    1380.0, 1390.0, 1400.0, 1410.0, 1420.0, 1430.0,
    1440.0, 1450.0, 1460.0, 1470.0, 1480.0, 1490.0,
    1500.0,
]
LEGACY_GAP = [
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.29153e-05,
    1.19734e-05, 1.11614e-05, 1.04004e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05,
    1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05,
    1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05,
    1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05,
    1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05,
    1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05,
    1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05,
    1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05,
    1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05, 1.00000e-05,
    1.03692e-05, 1.12328e-05, 1.24698e-05, 1.29808e-05, 1.37081e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05, 1.37572e-05,
    1.37572e-05,
]
LEGACY_PFC = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.7623, 5.7992, 9.4135,
    12.719, 15.904, 18.734, 18.898, 18.898, 18.898,
    18.898, 18.898, 18.898, 18.898, 18.898, 18.898,
    18.898, 18.898, 18.898, 18.898, 18.898, 18.898,
    18.898, 18.898, 18.898, 18.898, 18.898, 18.898,
    18.898, 18.898, 18.898, 18.898, 18.898, 18.898,
    18.898, 18.898, 18.898, 18.898, 18.898, 18.898,
    18.898, 18.898, 18.898, 18.898, 18.898, 18.898,
    18.898, 15.295, 11.701, 8.9789, 5.3113, 1.5958,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0,
]
LEGACY_RFO = [
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2180600e-03,
    4.2191600e-03, 4.2201100e-03, 4.2210100e-03, 4.2218900e-03, 4.2229500e-03, 4.2239100e-03,
    4.2247900e-03, 4.2256400e-03, 4.2263900e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03,
    4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03,
    4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03,
    4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03,
    4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03,
    4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03,
    4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03, 4.2264300e-03,
    4.2264300e-03, 4.2254700e-03, 4.2245200e-03, 4.2237900e-03, 4.2228200e-03, 4.2218400e-03,
    4.2210500e-03, 4.2200300e-03, 4.2185800e-03, 4.2179800e-03, 4.2171400e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03, 4.2170800e-03,
    4.2170800e-03,
]
LEGACY_RCI = [
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2309800e-03,
    4.2311400e-03, 4.2312800e-03, 4.2314100e-03, 4.2318900e-03, 4.2329500e-03, 4.2339100e-03,
    4.2347900e-03, 4.2356400e-03, 4.2363900e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03,
    4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03,
    4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03,
    4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03,
    4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03,
    4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03,
    4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03, 4.2364300e-03,
    4.2364300e-03, 4.2354700e-03, 4.2345200e-03, 4.2337900e-03, 4.2328200e-03, 4.2318400e-03,
    4.2314100e-03, 4.2312600e-03, 4.2310500e-03, 4.2309700e-03, 4.2308500e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03, 4.2308400e-03,
    4.2308400e-03,
]

# Optional: path to a legacy FRED Fortran output directory for live comparison.
# Set to a valid local path, or leave as None to use the embedded reference arrays below.
TRUE_LEGACY_DIR = None

_FORTRAN_NUM = re.compile(r'([+-]?\d+\.\d+)([+-]\d+)')


def _fix_fortran_float(s: str) -> str:
    return _FORTRAN_NUM.sub(r'\1E\2', s)


def _extract_single_value(line: str) -> float:
    values = [float(_fix_fortran_float(v)) for v in line[27:].split()]
    if not values:
        raise ValueError(f"Could not parse numeric value from line: {line.rstrip()}")
    return values[0]


def parse_legacy_outfrd(case_dir: str):
    times, gap, pfc, rfo, rci = [], [], [], [], []
    outfrd_files = sorted(glob.glob(os.path.join(case_dir, "outfrd*")))
    if not outfrd_files:
        raise FileNotFoundError(f"No outfrd files found in: {case_dir}")
    for fp in outfrd_files:
        vals = {}
        with open(fp, "r", encoding="utf-8") as fh:
            for line in fh:
                for key in ("time (s)", "gap (m)", "pfc (MPa)", "rfo (m)", "rci (m)"):
                    if line.startswith(key):
                        vals[key] = _extract_single_value(line)
        missing = [k for k in ("time (s)", "gap (m)", "pfc (MPa)", "rfo (m)", "rci (m)") if k not in vals]
        if missing:
            raise ValueError(f"Missing {', '.join(missing)} in {fp}")
        times.append(vals["time (s)"]); gap.append(vals["gap (m)"])
        pfc.append(vals["pfc (MPa)"]); rfo.append(vals["rfo (m)"]); rci.append(vals["rci (m)"])
    return (np.array(times), np.array(gap), np.array(pfc),
            np.array(rfo),   np.array(rci))


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
g = fred.FuelRodGeometry()
g.nf   = 3
g.nc   = 2
g.nz   = 1
g.rfi0 = 0.0
g.rfo0 = 4.2e-3
g.rci0 = 4.215e-3
g.rco0 = 4.9e-3
g.dz0  = [0.10]
g.vgp  = 1.0e-6
g.ruff = 5.0e-6
g.rufc = 5.0e-6
g.build()

# ---------------------------------------------------------------------------
# Materials / solver
# ---------------------------------------------------------------------------
fuel = fred.DummyFuelPellet()
clad = fred.DummyCladding()
gap  = fred.DummyGapMaterial()   # He gas bond — uses actual gap width + radiation

solver = fred.FredRodSolver(g, fuel, clad, gap)
solver.set_enable_heat_conduction(True)
solver.set_enable_stress_strain(True)

T0     = 700.0
GPRES0 = 0.1    # MPa, fill-gas pressure at T_REF = 293.15 K
solver.set_initial_temperature(T0)
solver.set_coolant_temperature([0.0, 1500.0], [T0, T0])
solver.set_power_density_history(
    [0.0, 100.0, 200.0, 600.0, 700.0, 1500.0],
    [0.0,   0.0, 5.0e8, 5.0e8,   0.0,    0.0],
)
solver.set_coolant_pressure(5.0)
solver.set_initial_gas_pressure(GPRES0)

solver.run(1500.0, 10.0)

# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
times = np.array(solver.time_points())
gap_w = np.array(solver.gap_width())
pfc   = np.array(solver.contact_pressure())
rfo   = np.array(solver.fuel_outer_radius())
rci   = np.array(solver.clad_inner_radius())

y_out = solver.y_out()
gpres = y_out[:, 0] if y_out.size > 0 else np.zeros_like(times)

thresh  = g.ruff + g.rufc
t_leg   = np.array(LEGACY_TIMES)
gap_leg = np.array(LEGACY_GAP)
pfc_leg = np.array(LEGACY_PFC)
rfo_leg = np.array(LEGACY_RFO)
rci_leg = np.array(LEGACY_RCI)

if TRUE_LEGACY_DIR and os.path.isdir(TRUE_LEGACY_DIR):
    try:
        t_leg, gap_leg, pfc_leg, rfo_leg, rci_leg = parse_legacy_outfrd(TRUE_LEGACY_DIR)
        print(f"Using true legacy outfrd data from: {TRUE_LEGACY_DIR}")
    except Exception as exc:
        print(f"Warning: failed to parse legacy data ({exc}); using embedded arrays.")
else:
    print("Using embedded legacy arrays (true legacy directory not found).")


def max_diff(t_new, v_new, t_ref, v_ref):
    return np.max(np.abs(np.interp(t_ref, t_new, v_new) - v_ref))


# Phase checks
phase1_open   = (gap_w[times <= 100.0] > thresh).all()
phase2_closed = (gap_w[(times >= 200.0) & (times <= 600.0)] <= thresh + 1e-9).all()
phase3_reopen = len(gap_w[times >= 800.0]) > 0 and (gap_w[times >= 800.0] > thresh).all()

print(f"\n=== FRED-ROD Gap Close/Reopen — Gas Bond vs Legacy FRED ===")
print(f"gpres0={GPRES0} MPa, T0={T0} K, pcool=5 MPa")
print(f"Gas pressure (DAE): {gpres.min():.4f} – {gpres.max():.4f} MPa")
print()
print(f"Phase 1 open-gap:   {'PASS' if phase1_open   else 'FAIL'}")
print(f"Phase 2 closed-gap: {'PASS' if phase2_closed else 'FAIL'}")
print(f"Phase 3 reopen:     {'PASS' if phase3_reopen else 'FAIL'}")
print()

mask_close  = (gap_w <= thresh) & (times > 100.0)
mask_reopen = (gap_w > thresh)  & (times > 600.0)
t_close  = times[mask_close].min()  if mask_close.any()  else float('nan')
t_reopen = times[mask_reopen].min() if mask_reopen.any() else float('nan')

leg_close  = t_leg[gap_leg <= thresh + 1e-10][0] if (gap_leg <= thresh + 1e-10).any() else float('nan')
leg_reopen = t_leg[(gap_leg > thresh) & (t_leg > 600.0)][0] if ((gap_leg > thresh) & (t_leg > 600.0)).any() else float('nan')

print(f"Gap closure:   FRED-ROD t = {t_close:.1f} s   |   Legacy t ≈ {leg_close:.0f} s  (Δ ≈ {t_close - leg_close:+.1f} s)")
if not np.isnan(t_reopen):
    print(f"Gap reopening: FRED-ROD t = {t_reopen:.1f} s  |   Legacy t ≈ {leg_reopen:.0f} s  (Δ ≈ {t_reopen - leg_reopen:+.1f} s)")
else:
    print(f"Gap reopening: FRED-ROD — did not reopen   |   Legacy t ≈ {leg_reopen:.0f} s")
print()

phases = [
    ("Phase 1 (0–100 s)",    t_leg <= 100),
    ("Phase 2 (200–600 s)",  (t_leg >= 200) & (t_leg <= 600)),
    ("Phase 3 (700–1500 s)", t_leg >= 700),
]
print(f"{'Region':<24}  {'Max|Δgap| µm':>14}  {'Max|Δpfc| MPa':>14}  {'Max|Δrfo| µm':>13}  {'Max|Δrci| µm':>13}")
print("-" * 85)
for label, mask in phases:
    if mask.sum() == 0:
        continue
    tl, gl, pl, rfl, rcl = t_leg[mask], gap_leg[mask], pfc_leg[mask], rfo_leg[mask], rci_leg[mask]
    print(f"{label:<24}  {max_diff(times, gap_w, tl, gl)*1e6:>14.4f}"
          f"  {max_diff(times, pfc, tl, pl):>14.4f}"
          f"  {max_diff(times, rfo, tl, rfl)*1e6:>13.4f}"
          f"  {max_diff(times, rci, tl, rcl)*1e6:>13.4f}")

print()
print("Overall max differences (whole run):")
print(f"  Max|Δgap| = {max_diff(times, gap_w, t_leg, gap_leg)*1e6:.4f} µm")
print(f"  Max|Δpfc| = {max_diff(times, pfc,   t_leg, pfc_leg):.4f} MPa")
print(f"  Max|Δrfo| = {max_diff(times, rfo,   t_leg, rfo_leg)*1e6:.4f} µm")
print(f"  Max|Δrci| = {max_diff(times, rci,   t_leg, rci_leg)*1e6:.4f} µm")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
plot_dir = os.path.join(os.path.dirname(__file__), 'plots')
os.makedirs(plot_dir, exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle(f'Gap Close/Reopen — Gas Bond (p₀={GPRES0} MPa) vs Legacy FRED', fontsize=12)

ax = axes[0, 0]
ax.plot(times, gap_w * 1e6, 'b-', lw=2, label='FRED-ROD')
ax.plot(t_leg, gap_leg * 1e6, 'r--', lw=1.5, label='Legacy FRED')
ax.axhline(thresh * 1e6, color='gray', ls=':', lw=1, label=f'threshold {thresh*1e6:.0f} µm')
ax.set_xlabel('Time (s)'); ax.set_ylabel('Gap width (µm)')
ax.set_title('Gap width'); ax.legend(); ax.grid(True)

ax = axes[0, 1]
ax.plot(times, pfc,    'b-',  lw=2,   label='FRED-ROD pfc')
ax.plot(t_leg, pfc_leg,'r--', lw=1.5, label='Legacy pfc')
ax.plot(times, gpres,  'g:',  lw=1.5, label='p_gas (DAE)')
ax.set_xlabel('Time (s)'); ax.set_ylabel('Pressure (MPa)')
ax.set_title('Contact pressure & gas pressure'); ax.legend(); ax.grid(True)

ax = axes[1, 0]
ax.plot(times, (rfo - rfo[0]) * 1e6, 'b-', lw=2, label='FRED-ROD')
ax.plot(t_leg, (rfo_leg - rfo_leg[0]) * 1e6, 'r--', lw=1.5, label='Legacy')
ax.set_xlabel('Time (s)'); ax.set_ylabel('Δrfo (µm)')
ax.set_title('Fuel outer radius change'); ax.legend(); ax.grid(True)

ax = axes[1, 1]
ax.plot(times, (rci - rci[0]) * 1e6, 'b-', lw=2, label='FRED-ROD')
ax.plot(t_leg, (rci_leg - rci_leg[0]) * 1e6, 'r--', lw=1.5, label='Legacy')
ax.set_xlabel('Time (s)'); ax.set_ylabel('Δrci (µm)')
ax.set_title('Clad inner radius change'); ax.legend(); ax.grid(True)

fig.tight_layout()
fig.savefig(os.path.join(plot_dir, 'results.png'), dpi=150)
plt.close(fig)
print(f"\nSaved {plot_dir}/results.png")
