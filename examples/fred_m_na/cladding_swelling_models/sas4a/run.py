#!/usr/bin/env python3
"""
Cladding void swelling — SAS4A model (Cswel.for icswel=2, added by
Timpano/Daniele 2024-08-06).

Runs the Timpano IF-AVG benchmark case (2176 d) with cladding void swelling
enabled and set to the SAS4A model: a piecewise-linear strain RATE, gated
on cladding dose >= 100 dpa, evaluated once at the cladding mid-wall
temperature and integrated (explicit Euler), applied uniformly across the
clad wall.

Outputs: sas4a.h5 (swelling/ecs, swelling/neuflue2, swelling/dose)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'build'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import _fred_m_na as fred
from common import run_case

run_case(fred, 'sas4a', os.path.dirname(os.path.abspath(__file__)))
