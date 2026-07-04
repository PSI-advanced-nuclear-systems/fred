"""
FRED-ROD advanced prototyping — Step 1: Python implementation.

Demonstrates the fission gas release gap conductance model implemented as a
pure-Python GapMaterial subclass.  No C++ recompilation is required.

Physics:
  Fission gas (Xe, Kr, He) released during irradiation mixes with the
  helium fill gas.  Xe and Kr have much lower thermal conductivity than He,
  so gap conductance decreases significantly with burnup.

  The model here uses the Waltar-Reynolds restructured-zone FGR fraction
  and a geometric-mean conductivity rule over the He/Kr/Xe mixture.

This file:
  1. Defines HeKrXePy — a Python GapMaterial with the burnup-dependent model.
  2. Prints the gas composition and gap conductance at several burnup levels.
  3. Runs three FRED-ROD simulations (bup = 0, 1, 5 at%) and compares the
     steady-state peak fuel temperature to show the thermal effect.

Once satisfied with the model, proceed to run_cpp.py which uses the compiled
C++ equivalent (fred.HeKrXe) with identical physics and no Python overhead.
"""

import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'build'))

import fred_rod as fred

# ===========================================================================
# Python GapMaterial: burnup-dependent He/Kr/Xe mixture
# ===========================================================================

class HeKrXePy(fred.GapMaterial):
    """
    Gap material whose conductivity reflects the He/Kr/Xe composition
    resulting from fission gas release at a given fuel burnup.

    Physics summary
    ---------------
    Step 1 — He fill-gas inventory at fabrication (ideal gas law):
        n_He0 = P_fill * V_free / (R * T_ref)
        Typical: P_fill=0.1 MPa, V_free=1.5 cm³, T_ref=293.15 K → n_He0 ≈ 6.15e-5 mol

    Step 2 — Fission gas generated and released (Waltar-Reynolds upper bound):
        fggen = 0.25 * (bup_atpct/100) * m_fuel / M_HM
                (0.25 gas atoms/fission, 10 g fuel, M_HM=238 g/mol)
        bup_MWd = bup_atpct * 9.4       (1 at% ≈ 9.4 MWd/kgHM, oxide)
        a   = 4.7 / bup_MWd * (1 - exp(-bup_MWd / 5.9))
        fgr = max(0,  1 - a)            (FGR fraction, restructured zone)
        fgrel = fgr * fggen

    Step 3 — Species moles (FRED.f90 gaphtc nuclear-physics composition):
        n_Xe   = 0.8846 * fgrel
        n_Kr   = 0.0769 * fgrel
        n_He_FG= 0.0385 * fgrel

    Step 4 — Mole fractions and geometric-mean conductivity rule:
        y_i    = n_i / n_total
        k_mix  = k_He^y_He * k_Kr^y_Kr * k_Xe^y_Xe
    """

    # Conductivity power-law fits  k(T) = A * T^B  [W/(m·K)]  (FRED.f90 gaphtc)
    _A_He, _B_He = 2.639e-3, 0.7085
    _A_Kr, _B_Kr = 8.247e-5, 0.8363
    _A_Xe, _B_Xe = 4.351e-5, 0.8618

    def __init__(self, bup_atpct: float = 0.0):
        super().__init__()
        bup_MWd = bup_atpct * 9.4

        # Waltar-Reynolds FGR fraction
        fgr = 0.0
        if bup_MWd > 0.0:
            a   = 4.7 / bup_MWd * (1.0 - math.exp(-bup_MWd / 5.9))
            fgr = max(0.0, 1.0 - a)

        # Molar balance
        R, T_ref        = 8.314, 293.15
        n_He0           = 0.1e6 * 1.5e-6 / (R * T_ref)     # fill-gas moles
        fggen           = 0.25 * (bup_atpct / 100.0) * 0.010 / 0.238
        fgrel           = fgr * fggen

        n_Xe            = 0.8846 * fgrel
        n_Kr            = 0.0769 * fgrel
        n_He_FG         = 0.0385 * fgrel
        n_total         = n_He0 + n_He_FG + n_Kr + n_Xe

        self.y_He = (n_He0 + n_He_FG) / n_total
        self.y_Kr = n_Kr / n_total
        self.y_Xe = n_Xe / n_total
        self.fgr  = fgr

    def gapConductivity(self, T: float) -> float:
        """Geometric-mean conductivity of the He/Kr/Xe mixture [W/(m·K)]."""
        k_He = self._A_He * T ** self._B_He
        k_Kr = self._A_Kr * T ** self._B_Kr
        k_Xe = self._A_Xe * T ** self._B_Xe
        return k_He ** self.y_He * k_Kr ** self.y_Kr * k_Xe ** self.y_Xe


# ===========================================================================
# Inspect composition and conductance at key burnup levels
# ===========================================================================

BURNUPS = [0.0, 0.5, 1.0, 2.0, 5.0]   # at%
T_GAP   = 700.0                         # K — representative gap temperature

print("=" * 72)
print("Gap gas composition and conductance vs burnup  (T_gap = 700 K)")
print("=" * 72)
print(f"  {'bup [at%]':>10}  {'FGR [%]':>8}  {'y_He':>6}  {'y_Kr':>6}  {'y_Xe':>6}"
      f"  {'k_gap [W/mK]':>13}  {'h_gap@50µm [W/m²K]':>20}")
print("  " + "-" * 75)
for bup in BURNUPS:
    g = HeKrXePy(bup)
    k = g.gapConductivity(T_GAP)
    h = k / 50.0e-6   # illustrative: k / 50 µm reference gap
    print(f"  {bup:>10.1f}  {g.fgr*100:>8.1f}  {g.y_He:>6.3f}  {g.y_Kr:>6.3f}"
          f"  {g.y_Xe:>6.3f}  {k:>13.4f}  {h:>20.1f}")

# ===========================================================================
# Geometry — single axial layer, solid oxide pellet
# ===========================================================================

def make_geometry():
    geom = fred.FuelRodGeometry()
    geom.nf   = 4
    geom.nc   = 3
    geom.nz   = 1
    geom.rfi0 = 0.0
    geom.rfo0 = 4.2e-3    # m
    geom.rci0 = 4.3e-3    # m  (100 µm as-fabricated gap)
    geom.rco0 = 4.9e-3    # m
    geom.dz0  = [0.10]    # m
    geom.vgp  = 1.5e-6    # m³
    geom.ruff = 5.0e-6    # m
    geom.rufc = 5.0e-6    # m
    geom.build()
    return geom

# ===========================================================================
# Solver factory
# ===========================================================================

T_COOLANT = 620.0    # K
POWER     = 3.0e8    # W/m³  (representative fast reactor linear power)

def make_solver(geom, gap_mat):
    s = fred.FredRodSolver(geom, fred.UO2(), fred.AIM1(), gap_mat)
    s.set_enable_heat_conduction(True)
    s.set_enable_stress_strain(True)
    s.set_initial_temperature(T_COOLANT)
    s.set_coolant_temperature([0.0, 2000.0], [T_COOLANT, T_COOLANT])
    # Slow power ramp to help IDA initialisation
    s.set_power_density_history([0.0, 50.0, 2000.0], [0.0, POWER, POWER])
    s.set_coolant_pressure(0.5)
    s.set_initial_gas_pressure(0.1)
    s.set_tolerances(1e-6, 1e-8)
    return s

# ===========================================================================
# Compare three burnup states
# ===========================================================================

CASES = [
    (0.0, "as-fabricated (pure He)"),
    (1.0, "moderate burnup (1 at%)"),
    (5.0, "high burnup    (5 at%)"),
]

geom    = make_geometry()
tend    = 1500.0
dtout   = 150.0

print()
print("=" * 72)
print(f"Steady-state fuel temperature vs burnup  (Q = {POWER/1e6:.0f} MW/m³, "
      f"T_cool = {T_COOLANT:.0f} K)")
print("=" * 72)
print(f"  {'Burnup':>22}  {'T_peak [K]':>11}  {'gap [µm]':>10}  "
      f"{'sigh_clad [MPa]':>16}")
print("  " + "-" * 66)

results = {}
for bup, label in CASES:
    gap = HeKrXePy(bup)
    s   = make_solver(geom, gap)
    s.run(tend, dtout)
    T_pk  = s.peak_fuel_temperature()[-1]
    gw    = s.gap_width()[-1]
    sigh  = s.clad_outer_hoop_stress()[-1]
    results[bup] = (T_pk, gw, sigh)
    print(f"  {label:>22}  {T_pk:>11.1f}  {gw*1e6:>10.2f}  {sigh:>16.3f}")

dT_1 = results[1.0][0] - results[0.0][0]
dT_5 = results[5.0][0] - results[0.0][0]
print(f"\n  Delta T_peak  (1 at% vs 0): {dT_1:+.1f} K")
print(f"  Delta T_peak  (5 at% vs 0): {dT_5:+.1f} K")

print("""
Notes:
  Higher burnup → more Xe/Kr in gap → lower k_gap → higher fuel temperature.
  The Waltar-Reynolds formula gives an upper-bound FGR (restructured-zone
  only); actual FGR at a given burnup will typically be lower.

  To compile this model into C++ and expose it to Python without overhead,
  see HeKrXe.hpp and run_cpp.py.
""")
