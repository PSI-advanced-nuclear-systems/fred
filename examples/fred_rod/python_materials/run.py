"""
FRED-ROD: Python-defined custom material properties — rapid-prototyping example.

Demonstrates how to subclass fred.FuelPelletMaterial and fred.CladdingMaterial
in Python to supply custom correlations without recompiling any C++.

Custom materials defined here:
  GenericOxideFuel   — simplified oxide fuel with analytic k(T) and cp(T)
  GenericSteelClad   — simplified ferritic steel with linear k(T) and E(T)

Scenario: heat conduction + stress-strain, ramp power 0→1.5e8 W/m3 over 50 s,
          hold for 950 s.  Coolant at 620 K.

After the run the peak fuel temperature and outer cladding hoop stress are
compared between the two custom materials and the built-in DummyFuelPellet /
DummyCladding for a quick sanity check.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'build'))

import fred_rod as fred

# ===========================================================================
# Custom fuel pellet — pure Python subclass of fred.FuelPelletMaterial
# ===========================================================================

T_REF = 293.15  # K — reference temperature for thermal expansion

class GenericOxideFuel(fred.FuelPelletMaterial):
    """
    Simplified oxide fuel material defined entirely in Python.

    Thermal conductivity:
        k(T) = 1 / (A + B*T) + C*T^3      [W/(m·K)]
    where A=0.045, B=2.8e-4, C=5.6e-11 gives ~6.5 W/mK at 600 K,
    decreasing to ~3 W/mK at 1500 K (qualitatively UO2-like).

    Heat capacity:
        cp(T) = 260 + 0.10*T               [J/(kg·K)]
    (linear fit, ~320 J/kgK at 600 K, increasing to ~410 J/kgK at 1500 K).

    Thermal expansion (linear):
        ε_th(T) = 1.0e-5 * (T - T_ref)    [-]

    Young's modulus (softens near melting):
        E(T) = E0 * (1 - 1.5e-4 * (T - T_ref))   [MPa]

    All remaining properties are constants typical for a dense oxide ceramic.
    """

    _A  = 0.045      # [m·K/W]  — phonon scattering coefficient
    _B  = 2.8e-4     # [m/W]    — phonon-phonon term
    _C  = 5.6e-11    # [W/(m·K^4)] — small-polaron / radiation term
    _E0 = 2.0e5      # [MPa]    — Young's modulus at T_ref (200 GPa)
    _rho_ref = 10500.0   # [kg/m³]  — as-fabricated density (~95.8 % TD)
    _rho_th  = 10960.0   # [kg/m³]  — theoretical density
    _T_melt  = 3100.0    # [K]

    def thermalConductivity(self, T):
        """k(T) = 1/(A + B*T) + C*T^3  [W/(m·K)]"""
        return 1.0 / (self._A + self._B * T) + self._C * T**3

    def heatCapacity(self, T):
        """cp(T) = 260 + 0.10*T  [J/(kg·K)]"""
        return 260.0 + 0.10 * T

    def thermalExpansionStrain(self, T):
        """ε_th(T) = 1e-5 * (T - T_ref)  [-]"""
        return 1.0e-5 * (T - T_REF)

    def youngsModulus(self, T, density):
        """E(T) = E0 * (1 - 1.5e-4*(T - T_ref))  [MPa]; density unused."""
        return self._E0 * max(0.0, 1.0 - 1.5e-4 * (T - T_REF))

    def poissonRatio(self):
        return 0.30

    def referenceDensity(self):
        return self._rho_ref

    def theoreticalDensity(self):
        return self._rho_th

    def meltingTemperature(self):
        return self._T_melt


# ===========================================================================
# Custom cladding — pure Python subclass of fred.CladdingMaterial
# ===========================================================================

class GenericSteelClad(fred.CladdingMaterial):
    """
    Simplified ferritic-martensitic steel cladding defined in Python.

    Thermal conductivity:
        k(T) = 23.5 + 0.005*(T - 273.15)   [W/(m·K)]
    (linear, typical range 23.5 at 0°C to ~31 at 1500°C).

    Heat capacity:
        cp(T) = 480 + 0.12*T                [J/(kg·K)]

    Thermal expansion (linear CTE = 11.5e-6 /K):
        ε_th(T) = 11.5e-6 * (T - T_ref)    [-]

    Young's modulus (linear softening):
        E(T) = 2.0e5 - 60.0*(T - T_ref)    [MPa]  (200 GPa → ~120 GPa at 1500 K)

    Meyer hardness (power-law):
        H_M(T) = 5.0e9 * T^(-0.2)          [Pa]
    """

    _k0  = 23.5     # [W/(m·K)]
    _dk  = 0.005    # [W/(m·K·°C)]
    _cp0 = 480.0    # [J/(kg·K)]
    _dcp = 0.12     # [J/(kg·K²)]
    _CTE = 11.5e-6  # [1/K]
    _E0  = 2.0e5    # [MPa]
    _dE  = 60.0     # [MPa/K]
    _rho_ref = 7800.0  # [kg/m³]

    def thermalConductivity(self, T):
        return self._k0 + self._dk * (T - 273.15)

    def heatCapacity(self, T):
        return self._cp0 + self._dcp * T

    def thermalExpansionStrain(self, T):
        return self._CTE * (T - T_REF)

    def youngsModulus(self, T):
        return max(1.0e4, self._E0 - self._dE * (T - T_REF))

    def poissonRatio(self):
        return 0.30

    def meyerHardness(self, T):
        """H_M(T) = 5e9 * T^(-0.2)  [Pa]"""
        return 5.0e9 * T**(-0.2)

    def referenceDensity(self):
        return self._rho_ref


# ===========================================================================
# Helper — print k(T) and cp(T) at representative temperatures
# ===========================================================================

def print_material_summary(fuel, clad, label):
    print(f"\n  {label}")
    print(f"  {'T [K]':>7}  {'k_fuel [W/mK]':>14}  {'cp_fuel [J/kgK]':>16}"
          f"  {'k_clad [W/mK]':>14}  {'E_clad [GPa]':>13}")
    print("  " + "-" * 70)
    for T in (400, 700, 1000, 1400):
        print(f"  {T:>7}  {fuel.thermal_conductivity(T):>14.3f}"
              f"  {fuel.heat_capacity(T):>16.1f}"
              f"  {clad.thermal_conductivity(T):>14.3f}"
              f"  {clad.youngs_modulus(T)/1e3:>13.1f}")


# ===========================================================================
# Geometry — shared for all runs
# ===========================================================================

def make_geometry():
    g = fred.FuelRodGeometry()
    g.nf   = 4
    g.nc   = 3
    g.nz   = 1
    g.rfi0 = 0.0
    g.rfo0 = 4.2e-3      # m
    g.rci0 = 4.3e-3      # m  (100 µm as-fab gap)
    g.rco0 = 4.9e-3      # m
    g.dz0  = [0.10]
    g.vgp  = 1.0e-6
    g.ruff = 5.0e-6
    g.rufc = 5.0e-6
    g.build()
    return g


# ===========================================================================
# Solver factory
# ===========================================================================

T0_K = 620.0  # K  initial / coolant temperature

def make_solver(geom, fuel, clad, gap):
    s = fred.FredRodSolver(geom, fuel, clad, gap)
    s.set_enable_heat_conduction(True)
    s.set_enable_stress_strain(True)
    s.set_initial_temperature(T0_K)
    s.set_coolant_temperature([0.0, 1000.0], [T0_K, T0_K])
    # Ramp power 0 → 1.5e8 W/m3 over 50 s, then hold — avoids IDACalcIC failure
    s.set_power_density_history([0.0, 50.0, 1000.0], [0.0, 1.5e8, 1.5e8])
    s.set_coolant_pressure(0.1)
    s.set_initial_gas_pressure(0.1)
    s.set_tolerances(1e-6, 1e-8)
    return s


# ===========================================================================
# Runs
# ===========================================================================

tend  = 1000.0
dtout = 100.0

g   = make_geometry()
gap = fred.DummyGapMaterial()

# -- Run A: Python custom materials --
print("=" * 60)
print("Run A: GenericOxideFuel + GenericSteelClad (Python materials)")
print("=" * 60)

fuel_py = GenericOxideFuel()
clad_py = GenericSteelClad()

print_material_summary(fuel_py, clad_py, "Property table at key temperatures")

s_py = make_solver(g, fuel_py, clad_py, gap)
s_py.run(tend, dtout)

times_py  = s_py.time_points()
T_peak_py = s_py.peak_fuel_temperature()
sigh_py   = s_py.clad_outer_hoop_stress()
gap_py    = s_py.gap_width()

print(f"\n  t_final = {times_py[-1]:.0f} s")
print(f"  Peak fuel temperature : {T_peak_py[-1]:.1f} K")
print(f"  Outer clad hoop stress: {sigh_py[-1]:.3f} MPa")
print(f"  Gap width             : {gap_py[-1]*1e6:.2f} µm")

# -- Run B: built-in Dummy materials for comparison --
print()
print("=" * 60)
print("Run B: DummyFuelPellet + DummyCladding (built-in C++ materials)")
print("=" * 60)

fuel_dummy = fred.DummyFuelPellet()
clad_dummy = fred.DummyCladding()

print_material_summary(fuel_dummy, clad_dummy, "Property table at key temperatures")

s_dummy = make_solver(g, fuel_dummy, clad_dummy, gap)
s_dummy.run(tend, dtout)

times_d   = s_dummy.time_points()
T_peak_d  = s_dummy.peak_fuel_temperature()
sigh_d    = s_dummy.clad_outer_hoop_stress()
gap_d     = s_dummy.gap_width()

print(f"\n  t_final = {times_d[-1]:.0f} s")
print(f"  Peak fuel temperature : {T_peak_d[-1]:.1f} K")
print(f"  Outer clad hoop stress: {sigh_d[-1]:.3f} MPa")
print(f"  Gap width             : {gap_d[-1]*1e6:.2f} µm")

# ===========================================================================
# Side-by-side comparison
# ===========================================================================

print()
print("=" * 60)
print("Comparison at steady state (t = {:.0f} s)".format(times_py[-1]))
print("=" * 60)
print(f"  {'Quantity':<35}  {'Python mat.':>12}  {'Dummy mat.':>12}")
print("  " + "-" * 65)
print(f"  {'Peak fuel temperature [K]':<35}  {T_peak_py[-1]:>12.1f}  {T_peak_d[-1]:>12.1f}")
print(f"  {'Outer clad hoop stress [MPa]':<35}  {sigh_py[-1]:>12.3f}  {sigh_d[-1]:>12.3f}")
print(f"  {'Gap width [µm]':<35}  {gap_py[-1]*1e6:>12.2f}  {gap_d[-1]*1e6:>12.2f}")

print("""
Notes:
  The GenericOxideFuel has lower thermal conductivity than the Dummy material
  (DummyFuelPellet uses a constant k=10 W/mK), so the peak fuel temperature
  is substantially higher despite the same power density.  The gap is open in
  both cases (as-fabricated gap = 100 µm).

  To define your own material: subclass fred.FuelPelletMaterial or
  fred.CladdingMaterial and override all abstract methods.  The instance is
  passed directly to FredRodSolver — no C++ recompilation required.
""")
