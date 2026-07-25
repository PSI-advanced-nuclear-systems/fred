#!/usr/bin/env python3
"""
Cladding wastage — PrecipitationKinetics model (Clanth.for port).

Layer-lumped lanthanide inventory + sqrt(D/t) precipitation-front law with
SAS-refitted coefficients.  Wastage grows only while the gap is in soft or
clos contact; the threshold-demo plot picks the layer with the latest
contact onset to show wastage switching on exactly at first soft contact.

Outputs: fred_m_na_pk.h5 (incl. per-axial-layer burnup/xwast_layer),
         wastage_threshold_demo.png
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'build'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import _fred_m_na as fred

from common import make_solver, plot_threshold_demo, TEND, DTOUT

THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
OUTPUT_H5 = os.path.join(THIS_DIR, "fred_m_na_pk.h5")

solver = make_solver(fred)
solver.set_cladding_wastage_model(fred.CladWastageModel.PrecipitationKinetics)
solver.set_output_file(OUTPUT_H5)
solver.run(TEND, DTOUT)
print("Run complete.")

ok = plot_threshold_demo(OUTPUT_H5,
                         os.path.join(THIS_DIR, 'wastage_threshold_demo.png'),
                         'PrecipitationKinetics: wastage gated on soft/clos contact')
sys.exit(0 if ok else 1)
