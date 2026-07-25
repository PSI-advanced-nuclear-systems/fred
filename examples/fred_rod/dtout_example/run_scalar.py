#!/usr/bin/env python3
"""
dtout as a SCALAR — constant output interval (the classic behaviour).

The rod and boundary conditions are the coupled thermo-mechanical case of
examples/fred_rod/heat_conduction_stress_strain, run to 200 s with
dtout = 10 s.  Postprocessing validates that the solver produced output at
exactly the requested times: t = 0, 10, 20, ..., 200 s.

Run:  python run_scalar.py
"""

import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'build'))
import fred_rod as fred

from common import make_solver, validate_output_times

TEND  = 200.0   # s
DTOUT = 10.0    # s — constant output interval

solver = make_solver(fred)
solver.run(tend=TEND, dtout=DTOUT)

# ── Validation ───────────────────────────────────────────────────────────────
times    = np.asarray(solver.time_points())
expected = np.arange(0.0, TEND + 0.5 * DTOUT, DTOUT)

ok = validate_output_times(times, expected, label=f"scalar dtout = {DTOUT} s")
sys.exit(0 if ok else 1)
