#!/usr/bin/env python3
"""
dtout as a VECTOR — variable output-interval schedule.

Motivation: the first seconds of a coupled thermo-mechanical run are the
stiffest part (power step, thermal expansion onset).  Forcing very small
output intervals early on also bounds SUNDIALS IDA's internal step there
(the time loop sets IDASetStopTime at every output time), which can help
the solver through numerical instabilities; once the transient smooths
out the intervals grow and the run costs no more than a constant-dtout one.

Schedule used here (sums to 205 s >= tend = 200 s):
    10 x 0.1 s   — very fine while the power step hits
     8 x 0.5 s   — fine during the early heat-up
    10 x 2.0 s   — medium while approaching steady state
     7 x 25.0 s  — coarse once nothing changes any more

The Python wrapper enforces np.sum(dtout) >= tend; if the schedule ran out
anyway, the C++ loop would repeat the last interval.

Postprocessing validates the produced output times against np.cumsum of the
requested schedule and writes dtout_validation.png comparing scalar vs
vector spacing.

Run:  python run_vector.py
"""

import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'build'))
import fred_rod as fred

from common import make_solver, validate_output_times

TEND = 200.0   # s
DTOUT_SCHEDULE = [0.1] * 10 + [0.5] * 8 + [2.0] * 10 + [25.0] * 7

assert np.sum(DTOUT_SCHEDULE) >= TEND   # wrapper checks this too

solver = make_solver(fred)
solver.run(tend=TEND, dtout=DTOUT_SCHEDULE)

# ── Validation ───────────────────────────────────────────────────────────────
times = np.asarray(solver.time_points())

# Expected output times: cumulative sum of the schedule, clipped at tend.
expected = np.concatenate([[0.0], np.cumsum(DTOUT_SCHEDULE)])
expected = expected[expected <= TEND + 1e-9]

ok = validate_output_times(times, expected, label="vector dtout schedule")

# ── Sanity check that the wrapper rejects an under-covering schedule ─────────
try:
    solver_bad = make_solver(fred)
    solver_bad.run(tend=TEND, dtout=[1.0, 2.0, 3.0])   # sums to 6 s << 200 s
    print("FAIL: under-covering schedule was not rejected")
    ok = False
except ValueError as e:
    print(f"OK:   under-covering schedule rejected as expected\n      ({e})")

# ── Plot: scalar vs vector spacing ───────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    solver_s = make_solver(fred)
    solver_s.run(tend=TEND, dtout=10.0)
    t_s = np.asarray(solver_s.time_points())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(t_s[1:],   np.diff(t_s),   'o-',  label='scalar dtout = 10 s')
    ax1.plot(times[1:], np.diff(times), 's--', label='vector schedule')
    ax1.set_xlabel('time [s]'); ax1.set_ylabel('output interval [s]')
    ax1.set_yscale('log'); ax1.legend(); ax1.set_title('Output spacing')

    ax2.plot(t_s,   solver_s.peak_fuel_temperature(), 'o-',  label='scalar')
    ax2.plot(times, solver.peak_fuel_temperature(),   's--', label='vector')
    ax2.set_xlabel('time [s]'); ax2.set_ylabel('peak fuel T [K]')
    ax2.legend(); ax2.set_title('Same physics, different sampling')

    fig.tight_layout()
    out_png = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'dtout_validation.png')
    fig.savefig(out_png, dpi=140)
    print(f"Plot written: {out_png}")
except ImportError:
    print("matplotlib not available — skipping plot")

sys.exit(0 if ok else 1)
