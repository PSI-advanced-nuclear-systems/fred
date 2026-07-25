#!/usr/bin/env python3
"""
FRED-M-Na, dtout as a SCALAR — constant output interval (classic behaviour).

Simplified U-Pu-Zr/HT-9/sodium rod run for 10 days with dtout = 1 day.
Postprocessing validates that outputs land at exactly t = 0, 1, 2, ..., 10 d.

Run:  python run_scalar.py
"""

import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'build'))
import _fred_m_na as fred

from common import make_solver, validate_output_times, D2S

TEND  = 10 * D2S   # s
DTOUT = 1 * D2S    # s — constant output interval (1 day)

solver = make_solver(fred)
solver.run(TEND, DTOUT)

# ── Validation ───────────────────────────────────────────────────────────────
times    = np.asarray(solver.times())
expected = np.arange(0.0, TEND + 0.5 * DTOUT, DTOUT)

ok = validate_output_times(times, expected, label="FRED-M-Na scalar dtout = 1 d")
sys.exit(0 if ok else 1)
