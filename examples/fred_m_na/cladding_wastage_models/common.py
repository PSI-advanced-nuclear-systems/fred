"""
Shared setup and plotting helpers for the cladding-wastage model examples.

Rod, power, and boundary conditions follow the Timpano SAS4A benchmark
(examples/fred_m_na/timpano_SAS_benchmark/base_irradiation): ESFR-SIMPLE
U-Pu-Zr/HT-9 geometry, the benchmark's 24-layer time-dependent power table
(parsed from its ref_input.txt, card 300000), SigmoidBurnup conductivity,
Subbotin HTC, 2176-day base irradiation, 1-day steps with a 1-second first
output (hot start).  This matches the configuration behind Timpano's
FRED-MFUEL lanthanide-wastage comparison (~43 um axial max at 2176 d).
"""

import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

D2S  = 86400.0
NF, NC, NZ = 11, 5, 24
RFO  = 3.863e-3
RCI  = RFO + 6.38e-4
RCO  = 5.025e-3
TCLAD0 = RCO - RCI              # as-fabricated cladding wall [m]
TEND  = 2176 * D2S
DTOUT = [1.0] + [D2S] * 2176

REF_INPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', 'timpano_SAS_benchmark', 'base_irradiation',
                         'ref_input.txt')


def load_benchmark_power(path=REF_INPUT):
    """Parse the 300000 power cards: [time, qv_z1..qv_z24] rows [W/m3].

    Returns (times, qqv_per_layer) ready for
    set_power_density_history_per_layer, with constant extrapolation to 2e8 s
    (same convention as the benchmark run.py).
    """
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if parts and parts[0] == '300000':
                rows.append([float(x) for x in parts[1:]])
    rows = np.array(rows)
    assert rows.shape[1] == NZ + 1, f"expected {NZ} power columns"
    times = rows[:, 0].tolist()
    qqv   = rows[:, 1:].T.tolist()          # [nz][ntimes]
    times.append(2.0e8)
    for j in range(NZ):
        qqv[j].append(qqv[j][-1])
    return times, qqv


def make_solver(fred):
    g = fred.FuelRodGeometry()
    g.nf   = NF
    g.nc   = NC
    g.nz   = NZ
    g.rfi0 = 0.0
    g.rfo0 = RFO
    g.rci0 = RCI
    g.rco0 = RCO
    g.dz0  = [0.03125] * NZ
    g.vgp  = 7.6375e-5
    g.ruff = 1.0e-6
    g.rufc = 1.0e-6
    g.build()

    fuel = fred.UPuZr(pu_weight_frac=0.1325, zr_weight_frac=0.10,
                      reference_density=15800.0)
    clad = fred.HT9(reference_density=7634.5)

    solver = fred.FredMNaSolver(g, fuel, clad)
    solver.set_grsis_data_mode(fred.GrsisDataMode.FEAST)
    solver.set_sodium_mode(fred.SodiumMode.TDependent)
    # SigmoidBurnup + benchmark power table: matches the Timpano SAS4A
    # benchmark setup used for the FRED-MFUEL wastage comparison.
    solver.set_conductivity_model(fred.ConductivityModel.SigmoidBurnup)
    p_times, p_qqv = load_benchmark_power()
    solver.set_power_density_history_per_layer(p_times, p_qqv)
    solver.set_coolant_channel(
        dhyd          = 3.887e-3,
        xarea         = 3.301e-5,
        flowr         = 0.162,
        T_inlet_times = [0.0, 2.0e8],
        T_inlet_vals  = [633.0, 633.0],
        corr          = fred.HtcCorrelation.Subbotin,
    )
    solver.set_coolant_pressure(0.3856)
    solver.set_initial_temperature(633.0)
    solver.set_initial_gas_pressure(0.101)
    solver.set_step_size(D2S)
    solver.set_hot_start(True)
    return solver


def read_wastage(h5path):
    """Return times [s], burnup [at%], xwast_layer (n, NZ) [m], contact (n, NZ)."""
    with h5py.File(h5path, 'r') as f:
        t   = f["time"][:]
        bu  = f["burnup/bup_atpct"][:]
        xw  = f["burnup/xwast_layer"][:]
        con = f["thermal/contact_state"][:]
    return t, bu, xw, con


def pick_demo_layer(t, contact):
    """Layer with the latest first-soft time (clearest threshold demo)."""
    onset = np.full(contact.shape[1], np.inf)
    for j in range(contact.shape[1]):
        idx = np.nonzero(contact[:, j] >= 1.0)[0]
        if idx.size:
            onset[j] = t[idx[0]]
    j = int(np.argmax(np.where(np.isfinite(onset), onset, -np.inf)))
    return j, onset[j]


def plot_threshold_demo(h5path, out_png, title):
    """Wastage vs time at the demo layer, with soft/clos onset markers —
    demonstrates that wastage grows only once the gap is in soft/clos."""
    t, bu, xw, con = read_wastage(h5path)
    j, _ = pick_demo_layer(t, con)

    i_soft = np.nonzero(con[:, j] >= 1.0)[0]
    i_clos = np.nonzero(con[:, j] >= 2.0)[0]
    t_soft = t[i_soft[0]] / D2S if i_soft.size else None
    t_clos = t[i_clos[0]] / D2S if i_clos.size else None

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t / D2S, xw[:, j] * 1e6, '-', lw=1.5,
            label=f'wastage, layer {j+1}/{con.shape[1]}')
    if t_soft is not None:
        ax.axvline(t_soft, color='C1', ls='--', lw=1.0)
        ax.text(t_soft + 15, 0.05 * ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0.1,
                f'first soft ({t_soft:.0f} d)', rotation=90,
                va='bottom', fontsize=9, color='C1')
    if t_clos is not None:
        ax.axvline(t_clos, color='C3', ls='--', lw=1.0)
        ax.text(t_clos + 15, 0.4 * max(ax.get_ylim()[1], 0.1),
                f'first clos ({t_clos:.0f} d)', rotation=90,
                va='bottom', fontsize=9, color='C3')
    ax.set_xlabel('time [days]')
    ax.set_ylabel('cladding wastage [µm]')
    ax.set_title(title)
    ax.grid(alpha=0.3); ax.legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    print(f"Plot written: {out_png}")

    # Threshold check: wastage must be exactly zero before first soft contact
    if i_soft.size:
        pre = np.max(np.abs(xw[:i_soft[0], j]))
        ok = pre == 0.0
        print(f"Threshold check (layer {j+1}): max wastage before soft contact "
              f"= {pre:.3e} m -> {'PASS' if ok else 'FAIL'}")
        return ok
    return True
