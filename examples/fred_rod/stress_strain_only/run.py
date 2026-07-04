"""
Stress-strain only — FRED-ROD benchmark with DummyCladding.

Physics:  heat conduction OFF  (temperature pinned at T0 = 293.15 K)
          stress-strain   ON

Scenario: coolant pressure ramps from 5 MPa down to 0.1 MPa over 100 s,
          while internal gas pressure stays at 0.1 MPa.
          As external pressure decreases the cladding expands radially.

Thin-wall hoop stress estimate (inner pressure p_i, outer pressure p_o):
  σ_θ ≈ (p_i * r_i  -  p_o * r_o) / t_clad   [MPa]
  Initial: σ_θ ≈ (0.1*4.3e-3 - 5.0*4.9e-3) / 0.6e-3 ≈ -40 MPa  (compressive)
  Final:   σ_θ ≈ (0.1*4.3e-3 - 0.1*4.9e-3) / 0.6e-3 ≈  -1 MPa
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'build'))

import fred_rod as fred

# ── Geometry ─────────────────────────────────────────────────────────────────
g = fred.FuelRodGeometry()
g.nf   = 3
g.nc   = 2
g.nz   = 1
g.rfi0 = 0.0
g.rfo0 = 4.2e-3
g.rci0 = 4.3e-3
g.rco0 = 4.9e-3
g.dz0  = [0.10]
g.vgp  = 1.0e-6
g.ruff = 5.0e-6
g.rufc = 5.0e-6
g.build()

# ── Materials ─────────────────────────────────────────────────────────────────
fuel = fred.DummyFuelPellet()
clad = fred.DummyCladding()
gap  = fred.DummyGapMaterial()

# ── Solver ───────────────────────────────────────────────────────────────────
solver = fred.FredRodSolver(g, fuel, clad, gap)

solver.set_enable_heat_conduction(False)   # temperature frozen at T0
solver.set_enable_stress_strain(True)

T0 = 293.15   # K  — T_REF, zero thermal expansion
solver.set_initial_temperature(T0)

# Coolant pressure ramp: 5 MPa → 0.1 MPa over 100 s
t_p  = [0.0, 100.0]
p_co = [5.0,   0.1]   # MPa
solver.set_coolant_pressure_history(t_p, p_co)
solver.set_initial_gas_pressure(0.1)   # MPa  reference fill-gas pressure at T_REF=293.15 K
solver.set_tolerances(1e-6, 1e-8)

# ── Run ──────────────────────────────────────────────────────────────────────
tend  = 100.0
dtout = 10.0
solver.run(tend, dtout)

# ── Results ──────────────────────────────────────────────────────────────────
times      = solver.time_points()
sigh_outer = solver.clad_outer_hoop_stress()   # MPa, outer clad node
rco        = solver.clad_outer_radius()         # m

rco0 = g.rco0

print(f"\n{'t [s]':>8}  {'p_cool [MPa]':>14}  {'σ_θ_outer [MPa]':>17}  {'Δr_co [µm]':>12}")
print("-" * 58)

import math
for k, t in enumerate(times):
    # interpolate applied coolant pressure at this output time
    p_interp = p_co[0] + (p_co[1] - p_co[0]) * min(t / 100.0, 1.0)
    delta_rco = (rco[k] - rco0) * 1e6   # µm
    print(f"{t:8.1f}  {p_interp:14.2f}  {sigh_outer[k]:17.3f}  {delta_rco:12.4f}")

print(f"\nInitial outer clad hoop stress (analytic thin-wall):")
p_i, p_o = 0.1, 5.0
r_i, r_o = g.rci0, g.rco0
t_c = r_o - r_i
sig_tw = (p_i * r_i - p_o * r_o) / t_c
print(f"  σ_θ ≈ (p_i*r_i - p_o*r_o) / t = ({p_i}*{r_i:.4f} - {p_o}*{r_o:.4f}) / {t_c:.4f} = {sig_tw:.2f} MPa")

# Old FRED legacy values (stress_strain_only, tend=100s, dtout=10s)
# Source: examples/primitive_cases/dummy_materials/stress_strain_only/
sigh_outer_legacy = [-37.6481, -33.8949, -30.1413, -26.3874, -22.6331,
                     -18.8784, -15.1234, -11.3680,  -7.61226, -3.85614, -0.0996587]
rco_legacy        = [4.89854e-3, 4.89869e-3, 4.89883e-3, 4.89898e-3, 4.89913e-3,
                     4.89927e-3, 4.89942e-3, 4.89956e-3, 4.89971e-3, 4.89985e-3, 4.90000e-3]

print(f"\nComparison with legacy FRED:")
print(f"{'t [s]':>8}  {'σ_θ [MPa]':>10}  {'Legacy [MPa]':>13}  {'Δr_co [µm]':>11}  {'Legacy [µm]':>12}")
print("-" * 62)
for k, t in enumerate(times):
    delta_rco_new    = (rco[k] - rco0) * 1e6
    delta_rco_legacy = (rco_legacy[k] - rco0) * 1e6
    print(f"{t:8.1f}  {sigh_outer[k]:10.3f}  {sigh_outer_legacy[k]:13.3f}  {delta_rco_new:11.4f}  {delta_rco_legacy:12.4f}")
