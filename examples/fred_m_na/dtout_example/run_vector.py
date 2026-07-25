#!/usr/bin/env python3
"""
FRED-M-Na, dtout as a VECTOR — variable output-interval schedule.

Because FRED-M-Na's fixed backward-Euler integrator takes its internal step
from the output interval when no explicit set_step_size() is given, a vector
dtout directly controls the physics step: very small steps through the
startup, growing once the rod reaches thermal equilibrium.

Schedule (sums to exactly 10 days):
     5 x 0.1 d — fine while the rod heats up
     3 x 0.5 d — medium
     4 x 2.0 d — coarse at steady state

Postprocessing validates the produced output times against np.cumsum of the
requested schedule.

Run:  python run_vector.py
"""

import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'build'))
import _fred_m_na as fred

from common import make_solver, validate_output_times, D2S

TEND = 10 * D2S
DTOUT_SCHEDULE = [0.1 * D2S] * 5 + [0.5 * D2S] * 3 + [2.0 * D2S] * 4

assert np.sum(DTOUT_SCHEDULE) >= TEND

solver = make_solver(fred)
solver.run(TEND, DTOUT_SCHEDULE)

# ── Validation ───────────────────────────────────────────────────────────────
times = np.asarray(solver.times())

expected = np.concatenate([[0.0], np.cumsum(DTOUT_SCHEDULE)])
expected = expected[expected <= TEND + 1e-6]

ok = validate_output_times(times, expected,
                           label="FRED-M-Na vector dtout schedule")

# ── Plot: scalar vs vector spacing ───────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    solver_s = make_solver(fred)
    solver_s.run(TEND, 1 * D2S)
    t_s = np.asarray(solver_s.times())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(t_s[1:] / D2S,   np.diff(t_s) / D2S,   'o-',
             label='scalar dtout = 1 d')
    ax1.plot(times[1:] / D2S, np.diff(times) / D2S, 's--',
             label='vector schedule')
    ax1.set_xlabel('time [d]'); ax1.set_ylabel('output interval [d]')
    ax1.set_yscale('log'); ax1.legend(); ax1.set_title('Output spacing')

    ax2.plot(t_s / D2S,   solver_s.peak_fuel_temperature(), 'o-',
             label='scalar')
    ax2.plot(times / D2S, solver.peak_fuel_temperature(),   's--',
             label='vector')
    ax2.set_xlabel('time [d]'); ax2.set_ylabel('peak fuel T [K]')
    ax2.legend(); ax2.set_title('Same physics, different sampling')

    fig.tight_layout()
    out_png = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'dtout_validation.png')
    fig.savefig(out_png, dpi=140)
    print(f"Plot written: {out_png}")
except ImportError:
    print("matplotlib not available — skipping plot")

sys.exit(0 if ok else 1)
