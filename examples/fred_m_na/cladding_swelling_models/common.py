"""
Shared setup for the cladding void-swelling model examples (Hofmann vs
SAS4A, FredMNaCladSwelling.hpp).

Reuses the Timpano IF-AVG benchmark deck (geometry/coolant/power) parsed by
the sibling lanthanide_wastage/common.py, since the SAS4A model's fast-flux
calibration (fltpow, dconv) is itself an IF-AVG-specific Serpent fit —
using a different pin/power case without re-deriving fltpow would misstate
the neutron flux the swelling models see.
"""

import os
import sys
import importlib.util
import numpy as np
import h5py

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
LW_DIR   = os.path.join(THIS_DIR, '..', 'timpano_SAS_benchmark', 'lanthanide_wastage')

# Load the sibling lanthanide_wastage/common.py (IF-AVG deck parser) under a
# distinct module name to avoid colliding with this file's own "common".
_spec = importlib.util.spec_from_file_location(
    'lw_common', os.path.join(LW_DIR, 'common.py'))
_lw_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_lw_common)

D2S = _lw_common.D2S
DECK_DIR_OVERRIDE = os.path.join(LW_DIR, 'decks')


def load_if_avg_deck():
    _lw_common.DECK_DIR = DECK_DIR_OVERRIDE
    return _lw_common.parse_deck('if_avg')


def make_solver(fred, deck, swelling_model):
    g = fred.FuelRodGeometry()
    g.nf   = deck['nf']
    g.nc   = deck['nc']
    g.nz   = deck['nz']
    g.rfi0 = deck['rfi']
    g.rfo0 = deck['rfo']
    g.rci0 = deck['rfo'] + deck['dgap']
    g.rco0 = deck['rco']
    g.dz0  = [deck['dz0']] * deck['nz']
    g.vgp  = deck['vgp']
    g.ruff = deck['ruff']
    g.rufc = deck['rufc']
    g.build()

    fuel = fred.UPuZr(pu_weight_frac=deck['pu'], zr_weight_frac=deck['zr'],
                      reference_density=deck['fden'])
    clad = fred.HT9(reference_density=deck['cden'])

    solver = fred.FredMNaSolver(g, fuel, clad)
    solver.set_grsis_data_mode(fred.GrsisDataMode.FEAST)
    solver.set_sodium_mode(fred.SodiumMode.TDependent)
    solver.set_conductivity_model(fred.ConductivityModel.SigmoidBurnup)
    solver.set_power_density_history_per_layer(deck['power_times'],
                                               deck['power_qqv'])
    solver.set_coolant_channel(
        dhyd          = deck['dhyd'],
        xarea         = deck['xarea'],
        flowr         = deck['flowr'],
        T_inlet_times = [0.0, deck['tend'] + 1.0e7],
        T_inlet_vals  = [deck['tcoolin'], deck['tcoolin']],
        corr          = fred.HtcCorrelation.Subbotin,
    )
    solver.set_coolant_pressure(deck['press_Pa'] / 1.0e6)
    solver.set_initial_temperature(deck['tcoolin'])
    solver.set_initial_gas_pressure(deck['pin_Pa'] / 1.0e6)
    solver.set_step_size(D2S)
    solver.set_hot_start(True)

    # Cladding void swelling — off by default in the solver; this example
    # exists specifically to exercise it.
    solver.set_enable_clad_swelling(True)
    solver.set_clad_swelling_model(swelling_model)
    # fltpow/dconv default to the Timpano IF-AVG Serpent fit already, but
    # set explicitly here for clarity/traceability.
    solver.set_fast_flux_calibration(6.13209e14, 4.0)
    solver.set_clad_swelling_dose_threshold(100.0)
    return solver


def run_case(fred, model_name, this_dir, tend_days=2176):
    """model_name: 'hofmann' or 'sas4a'."""
    deck = load_if_avg_deck()
    model = (fred.CladSwellingModel.Hofmann if model_name == 'hofmann'
             else fred.CladSwellingModel.SAS4A)
    solver = make_solver(fred, deck, model)

    tend = min(deck['tend'], tend_days * D2S)
    n_days = int(np.ceil(tend / D2S))
    dtout = [1.0] + [D2S] * n_days

    h5 = os.path.join(this_dir, f'{model_name}.h5')
    solver.set_output_file(h5)
    solver.run(tend, dtout)
    print(f"[{model_name}] done -> {os.path.basename(h5)}")
    return h5


def read_swelling(h5path, nz, nc, layer=11, node=0):
    """layer/node are 0-based; legacy's OUTPUT_AXIAL_LAYER 12 -> layer=11,
    ecs(1,...) -> node=0 (innermost clad interval)."""
    with h5py.File(h5path, 'r') as f:
        t    = f['time'][:]
        ecs  = f['swelling/ecs'][:].reshape(-1, nz, nc)
        dose = f['swelling/dose'][:]
        neuflue2 = f['swelling/neuflue2'][:]
    return t, ecs[:, layer, node], dose[:, layer], neuflue2[:, layer]
