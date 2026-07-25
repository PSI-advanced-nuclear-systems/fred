#!/usr/bin/env python3
"""Timpano lanthanide-wastage benchmark, case of_avg: runs both wastage
models (PrecipitationKinetics + LaTracking) with deck-parsed conditions."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'build'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fred_m_na as fred
from common import run_case

run_case(fred, 'of_avg', os.path.dirname(os.path.abspath(__file__)))
