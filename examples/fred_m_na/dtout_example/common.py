"""
Shared rod setup and output-time validation for the FRED-M-Na dtout examples.

Rod and boundary conditions are a simplified version of
examples/fred_m_na/timpano_SAS_benchmark (U-Pu-Zr / HT-9 / sodium, ESFR-SIMPLE
geometry) with uniform constant power, run for a few days only.
"""

import numpy as np

D2S = 86400.0   # days -> seconds


def make_solver(fred):
    """Fresh solver per run (solver state is not resettable after run())."""
    NF, NC, NZ = 11, 5, 24
    RFO = 3.863e-3

    g = fred.FuelRodGeometry()
    g.nf   = NF
    g.nc   = NC
    g.nz   = NZ
    g.rfi0 = 0.0
    g.rfo0 = RFO
    g.rci0 = RFO + 6.38e-4
    g.rco0 = 5.025e-3
    g.dz0  = [0.03125] * NZ
    g.vgp  = 7.6375e-5
    g.ruff = 1.0e-6
    g.rufc = 1.0e-6
    g.build()

    fuel = fred.UPuZr(pu_weight_frac=0.1325, zr_weight_frac=0.10,
                      reference_density=15800.0)
    clad = fred.HT9(reference_density=7634.5)

    solver = fred.FredMNaSolver(g, fuel, clad)
    solver.set_grsis_data_mode(fred.GrsisDataMode.FEAST)
    solver.set_sodium_mode(fred.SodiumMode.TDependent)
    solver.set_conductivity_model(fred.ConductivityModel.SigmoidBurnup)
    solver.set_power_density_history([0.0, 2.0e8], [1.0e9, 1.0e9])
    solver.set_coolant_channel(
        dhyd          = 3.887e-3,
        xarea         = 3.301e-5,
        flowr         = 0.162,
        T_inlet_times = [0.0, 2.0e8],
        T_inlet_vals  = [633.0, 633.0],
        corr          = fred.HtcCorrelation.Subbotin,
    )
    solver.set_coolant_pressure(0.3856)
    solver.set_initial_temperature(633.0)
    solver.set_initial_gas_pressure(0.101)
    # NOTE: no set_step_size() — the internal backward-Euler step then follows
    # the output interval, so a vector dtout also refines the physics steps.
    return solver


def validate_output_times(times, expected, label, tol=1e-6):
    """Compare produced output times against the requested schedule.

    Checks count, absolute placement, and interval spacing (np.diff).
    Prints a PASS/FAIL report; returns True on PASS.
    """
    print(f"\n=== Validation: {label} ===")
    print(f"  requested output points : {len(expected)}")
    print(f"  produced  output points : {len(times)}")

    if len(times) != len(expected):
        print(f"  FAIL: count mismatch")
        print(f"        produced: {np.array2string(np.asarray(times), precision=3)}")
        print(f"        expected: {np.array2string(np.asarray(expected), precision=3)}")
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
