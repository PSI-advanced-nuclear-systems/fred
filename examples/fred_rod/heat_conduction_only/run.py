#!/usr/bin/env python3
"""
Heat conduction only — FRED-ROD benchmark with DummyFuelPellet and DummyCladding.

Physics:  heat conduction ON,  stress-strain OFF
Scenario: step power of 1e8 W/m3 applied at t=0; coolant held at 700 K.
          Integration runs to t=10 s (several thermal time constants).

Expected result: fuel centreline heats from 700 K to ~837 K at steady state.
  ΔT_fuel = Q * R_fuel² / (4 * k_fuel) = 1e8 * (4.2e-3)² / (4*10) ≈ 44 K  (fuel only)
  Additional ΔT across gap + cladding accounts for the remaining ~93 K.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'build'))

import fred_rod as fred

# ── Geometry ─────────────────────────────────────────────────────────────────
g = fred.FuelRodGeometry()
g.nf   = 3
g.nc   = 2
g.nz   = 1
g.rfi0 = 0.0        # solid pellet
g.rfo0 = 4.2e-3     # m
g.rci0 = 4.3e-3     # m   (100 µm gap)
g.rco0 = 4.9e-3     # m
g.dz0  = [0.10]     # m   (10 cm axial slice)
g.vgp  = 1.0e-6     # m3  (gas plenum)
g.ruff = 5.0e-6     # m   fuel surface roughness
g.rufc = 5.0e-6     # m   clad inner roughness
g.build()

# ── Materials ─────────────────────────────────────────────────────────────────
# Dummy material properties are defined in src/apps/fred_rod/{cladding,fuel,gap}material/Dummy{Cladding,FuelPellet,Gap}.cpp
# For Fred-Rod, new materials can be implemented two ways: (1) implement a new C++ class derived from fred::Material and recompile, or (2) use the Python interface to define a new material class derived from fred.Material.
fuel = fred.DummyFuelPellet() # constant properties k=10 W/mK, rho=10000 kg/m3, cp=100 J/kgK. 
clad = fred.DummyCladding()   # k=10 W/mK, cp=100 J/kgK
gap  = fred.DummyGapMaterial() # h=

# ── Solver ───────────────────────────────────────────────────────────────────
solver = fred.FredRodSolver(g, fuel, clad, gap)

# Physics toggles
solver.set_enable_heat_conduction(True)
solver.set_enable_stress_strain(False)

# Boundary conditions
T0 = 700.0   # K
solver.set_initial_temperature(T0)
solver.set_coolant_temperature([0.0, 1e6], [T0, T0])    # constant coolant BC
solver.set_power_density_history([0.0, 1e6], [1.0e8, 1.0e8])  # constant power

solver.set_coolant_pressure(5.0)   # MPa (not used by heat-only but required)
solver.set_initial_gas_pressure(0.1) # MPa
solver.set_tolerances(1e-6, 1e-8)

# ── Run ──────────────────────────────────────────────────────────────────────
tend  = 10.0   # s  (~30 × thermal time constant τ ≈ 0.3 s)
dtout = 0.5    # s
solver.run(tend, dtout)

# ── Results ──────────────────────────────────────────────────────────────────
times  = solver.time_points()
T_peak = solver.peak_fuel_temperature()
T_arr  = solver.temperatures()   # shape (n_times, nz*(nf+nc))

nf, nc, nz = g.nf, g.nc, g.nz
stride = nf + nc   # per axial layer (nz=1)

print(f"\n{'t [s]':>8}  {'T_center [K]':>14}  {'T_outer_fuel [K]':>18}  {'T_outer_clad [K]':>18}")
print("-" * 65)
for k, t in enumerate(times):
    T_row = T_arr[k]
    T_center     = T_row[0]          # fuel node 0 (innermost)
    T_outer_fuel = T_row[nf - 1]    # fuel node nf-1 (outer surface)
    T_outer_clad = T_row[nf + nc - 1]  # outer clad surface
    print(f"{t:8.2f}  {T_center:14.2f}  {T_outer_fuel:18.2f}  {T_outer_clad:18.2f}")

dT_ss = T_peak[-1] - T0
print(f"\nSteady-state ΔT (centreline above coolant): {dT_ss:.2f} K")

# print out all fuel and cladding temperatures at final time step
# verify that the final temperature matches legacy FRED results of : Tfuel: 842.531, 831.506, 798.431,  Tclad: 711.504, 700.  

Tfuel_legacy = [842.531, 831.506, 798.431]
Tclad_legacy = [711.504, 700.000]

print("\nFinal temperatures at t = {:.2f} s:".format(times[-1]))
print(f"{'Node':>6}  {'T [K]':>10}  {'Legacy FRED [K]':>12}  {'Material':>12}")
print("-" * 42)
for i in range(nf):
    T_fuel = T_arr[-1, i]
    print(f"{i:6d}  {T_fuel:10.2f} {Tfuel_legacy[i]:10.2f} {'Fuel':>12}")
for j in range(nc):
    T_clad = T_arr[-1, nf + j]
    print(f"{nf + j:6d}  {T_clad:10.2f} {Tclad_legacy[j]:10.2f} {'Cladding':>12}")
