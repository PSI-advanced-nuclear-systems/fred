"""
Shared rod setup and output-time validation for the dtout examples.

The rod and boundary conditions mirror
examples/fred_rod/heat_conduction_stress_strain: coupled thermo-elastic run
with a 50 s power ramp to 1e8 W/m3, coolant at 700 K, 5 MPa.
"""

import numpy as np


def make_solver(fred):
    """Fresh solver per run (solver state is not resettable after run())."""
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

    fuel = fred.DummyFuelPellet()
    clad = fred.DummyCladding()
    gap  = fred.DummyGapMaterial()

    solver = fred.FredRodSolver(g, fuel, clad, gap)
    solver.set_enable_heat_conduction(True)
    solver.set_enable_stress_strain(True)

    T0 = 700.0
    solver.set_initial_temperature(T0)
    solver.set_coolant_temperature_history([0.0, 1e6], [T0, T0])
    solver.set_power_density_history([0.0, 50.0, 100.0, 1e6],
                                     [0.0, 1.0e8, 1.0e8, 1.0e8])
    solver.set_coolant_pressure_history([0.0, 1e6], [5.0, 5.0])
    solver.set_initial_gas_pressure(0.1)
    solver.set_tolerances(1e-6, 1e-8)
    return solver


def validate_output_times(times, expected, label, tol=1e-9):
    """Compare produced output times against the requested schedule.

    Checks count, absolute placement, and interval spacing (np.diff).
    Prints a PASS/FAIL report; returns True on PASS.
    """
    print(f"\n=== Validation: {label} ===")
    print(f"  requested output points : {len(expected)}")
    print(f"  produced  output points : {len(times)}")

    if len(times) != len(expected):
        print(f"  FAIL: count mismatch")
        print(f"        produced: {np.array2string(times,    precision=3)}")
        print(f"        expected: {np.array2string(expected, precision=3)}")
        return False

    err_t  = np.max(np.abs(times - expected))
    err_dt = np.max(np.abs(np.diff(times) - np.diff(expected)))
    print(f"  max |t_produced - t_requested|      : {err_t:.3e} s")
    print(f"  max |dt_produced - dt_requested|    : {err_dt:.3e} s")

    if err_t < tol and err_dt < tol:
        print(f"  PASS: output times match the requested schedule")
        return True
    print(f"  FAIL: output times deviate from the requested schedule")
    return False
