#!/usr/bin/env python3
"""
Cladding void swelling — Hofmann (1985) model (Cswel.for icswel=1).

Runs the Timpano IF-AVG benchmark case (2176 d) with cladding void swelling
enabled and set to the Hofmann model: strain is a direct function of
cumulative fast-neutron fluence, evaluated per clad node at that node's own
temperature.

Outputs: hofmann.h5 (swelling/ecs, swelling/neuflue2, swelling/dose)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'build'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import _fred_m_na as fred
from common import run_case

run_case(fred, 'hofmann', os.path.dirname(os.path.abspath(__file__)))
