"""
Timpano lanthanide-wastage benchmark: shared deck parser and solver builder.

Three cases from timpano_fred/02.benchmark (decks copied into decks/):
  IF-AVG  : inner-fuel pin, core-average power   (0.75 m active, flowr 0.162)
  IF-PEAK : inner-fuel pin, peak power           (2-day startup ramp, 1.21e9 peak)
  OF-AVG  : outer-fuel pin, core-average power   (0.95 m active: dz0=0.0396,
            smaller plenum 6.36e-5 m3, flowr 0.120, ~60% power)

Everything geometry/coolant/power/tend is parsed from the legacy deck so the
platform runs stay faithful to what Timpano ran.  NOTE: the benchmark decks
specify FUEL_THK 14, a conductivity option that does not exist in the
available Flamb.for source; the platform runs use SigmoidBurnup (f=3, the
ESFR-SIMPLE fit) consistently for all three cases.
"""

import os
import numpy as np
import h5py

D2S = 86400.0
DECK_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'decks')

CASES = {
    'if_avg':  'if_avg.inp',
    'if_peak': 'if_peak.inp',
    'of_avg':  'of_avg.inp',
}


def parse_deck(case):
    """Parse geometry, coolant, time, and power cards from a legacy deck."""
    d = {}
    power_rows = []
    with open(os.path.join(DECK_DIR, CASES[case])) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == '000003':     # dz0, nz, nfr
                d['dz0'], d['nz'] = float(p[1]), int(p[2])
            elif p[0] == '100001':   # fmat fden pucont zrcont rfi rfo ruff stoch nf
                d['fden'], d['pu'], d['zr'] = float(p[2]), float(p[3]), float(p[4])
                d['rfi'], d['rfo'], d['ruff'] = float(p[5]), float(p[6]), float(p[7])
                d['nf'] = int(p[9])
            elif p[0] == '100002':   # gmat dgap pin vpl
                d['dgap'], d['pin_Pa'], d['vgp'] = float(p[2]), float(p[3]), float(p[4])
            elif p[0] == '100003':   # cmat rco cden rufc nc
                d['rco'], d['cden'], d['rufc'] = float(p[2]), float(p[3]), float(p[4])
                d['nc'] = int(p[5])
            elif p[0] == '500000':   # ctype dhyd xarea press flowr tcoolin
                d['dhyd'], d['xarea'] = float(p[2]), float(p[3])
                d['press_Pa'], d['flowr'], d['tcoolin'] = \
                    float(p[4]), float(p[5]), float(p[6])
            elif p[0] == '600000':   # tstart tend dtout
                d['tend'] = float(p[2])
            elif p[0] == '300000':   # time qv_1..qv_nz
                power_rows.append([float(x) for x in p[1:]])
    rows = np.array(power_rows)
    assert rows.shape[1] == d['nz'] + 1, f"{case}: power columns != nz"
    times = rows[:, 0].tolist()
    qqv   = rows[:, 1:].T.tolist()
    # constant extrapolation past the table end
    times.append(d['tend'] + 1.0e7)
    for j in range(d['nz']):
        qqv[j].append(qqv[j][-1])
    d['power_times'], d['power_qqv'] = times, qqv
    return d


def make_solver(fred, deck):
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
    return solver


def run_case(fred, case, this_dir, models=('pk', 'la')):
    """Run one case with the requested wastage models; returns H5 paths."""
    deck = parse_deck(case)
    n_days = int(np.ceil(deck['tend'] / D2S))
    dtout = [1.0] + [D2S] * n_days
    outs = {}
    for model_name in models:
        h5 = os.path.join(this_dir, f'{case}_{model_name}.h5')
        solver = make_solver(fred, deck)
        if model_name == 'pk':
            solver.set_cladding_wastage_model(
                fred.CladWastageModel.PrecipitationKinetics)
        else:
            solver.set_cladding_wastage_model(fred.CladWastageModel.LaTracking)
            # MFUEL's unconditional Dirichlet sink (transport through the
            # sodium bond): active from t=0, not gated on contact.
            solver.set_la_sink_requires_contact(False)
        solver.set_output_file(h5)
        solver.run(deck['tend'], dtout)
        outs[model_name] = h5
        print(f"[{case}/{model_name}] done -> {os.path.basename(h5)}")
    return outs


def axial_wastage_at(h5path, t_target_s):
    """Return (layer indices 1..nz, wastage [m] per layer) at ~t_target_s."""
    with h5py.File(h5path, 'r') as f:
        t  = f['time'][:]
        xw = f['burnup/xwast_layer'][:]
    i = int(np.argmin(np.abs(t - t_target_s)))
    return np.arange(1, xw.shape[1] + 1), xw[i], t[i]
