"""
FRED-M-Na — user-facing Python module.

U-Pu-Zr metallic fuel in sodium-bonded HT-9 cladding. Extends the FRED-ROD
thermal/mechanical solver with: Zr redistribution, cladding wastage, GRSIS
bubble swelling, fission gas, and a subchannel energy balance for coolant
temperature (Robin BC at cladding outer surface).

Typical usage
-------------
    import sys; sys.path.insert(0, '/path/to/fred_platform/build')
    import fred_m_na as fred

    g = fred.FuelRodGeometry()
    g.nf = 5; g.nc = 3; g.nz = 10
    g.rfi0 = 0.0;     g.rfo0 = 2.15e-3
    g.rci0 = 2.285e-3; g.rco0 = 2.67e-3
    g.dz0 = [0.085] * 10; g.vgp = 5e-7
    g.build()

    fuel   = fred.UPuZr(pu_wf=0.19, zr_wf=0.10)
    clad   = fred.HT9()
    solver = fred.FredMNaSolver(g, fuel, clad)

    solver.set_initial_temperature(700.0)
    solver.set_power_density_history([0, 1e8], [3e8, 3e8])
    solver.set_coolant_channel(
        dhyd=3.7e-3, xarea=1.5e-4, flowr=0.8,
        T_inlet_times=[0, 1e8], T_inlet_vals=[670.0, 670.0],
    )
    solver.set_initial_gas_pressure(0.1)
    solver.run(tend=1e8, dtout=1e7, output_file='results.h5')

    r = fred.read_results('results.h5')
"""

import numpy as np

import _fred_m_na as _cpp   # compiled C++ extension

from fred_geometry import FuelRodGeometry
from fred_materials import FuelPellet, Cladding, GapFill, _validate_density
from fred_solver_base import FredSolverBaseSubchannel


# ---------------------------------------------------------------------------
# MNA-specific material wrappers
# ---------------------------------------------------------------------------

class UPuZr(FuelPellet):
    """
    U-Pu-Zr metallic fuel pellet (Karahan 2009 correlations).

    Parameters
    ----------
    pu_wf             : Pu weight fraction [-] (e.g. 0.19 for 19 wt%)
    zr_wf             : Zr weight fraction [-] (e.g. 0.10 for 10 wt%)
    reference_density : as-fabricated density [kg/m³]; None → 75% of theoretical
    """
    def __init__(self, pu_wf, zr_wf, reference_density=None):
        if pu_wf <= 0 or pu_wf >= 1:
            raise ValueError(f"pu_wf must be in (0, 1), got {pu_wf!r}")
        if zr_wf <= 0 or zr_wf >= 1:
            raise ValueError(f"zr_wf must be in (0, 1), got {zr_wf!r}")
        if pu_wf + zr_wf >= 1:
            raise ValueError(
                f"pu_wf + zr_wf must be < 1 (remainder is U), "
                f"got {pu_wf} + {zr_wf} = {pu_wf + zr_wf}"
            )
        _validate_density("UPuZr", reference_density)
        self._pu = pu_wf
        self._zr = zr_wf
        self._ref_dens = reference_density
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = (cpp_module.UPuZr(self._pu, self._zr, self._ref_dens)
                             if self._ref_dens is not None
                             else cpp_module.UPuZr(self._pu, self._zr))
        return self._cpp_obj


class HT9(Cladding):
    """
    HT-9 ferritic-martensitic steel cladding (Hofmann 1985, Karahan 2007).

    Parameters
    ----------
    reference_density : as-fabricated density [kg/m³] (default 7750 kg/m³)
    """
    def __init__(self, reference_density=None):
        _validate_density("HT9", reference_density)
        self._ref_dens = reference_density
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = (cpp_module.HT9(self._ref_dens)
                             if self._ref_dens is not None else cpp_module.HT9())
        return self._cpp_obj


# ---------------------------------------------------------------------------
# FredMNaSolver
# ---------------------------------------------------------------------------

_CONDUCTIVITY_MODELS = {
    'DetailedNaSodium': None,
    'EmpiricalBurnup':  None,
    'EsfrSimple':       None,
}

_SODIUM_MODES = {'TDependent', 'Constant'}


class FredMNaSolver(FredSolverBaseSubchannel):
    """
    FRED-M-Na solver: metallic U-Pu-Zr fuel in sodium, HT-9 cladding.

    Coolant temperature is computed by a subchannel energy balance; the boundary
    condition at the cladding outer surface is Robin (heat-flux based).
    Use set_coolant_channel() — not set_coolant_temperature_history().
    """

    def __init__(self, geometry, fuel, cladding):
        if not isinstance(fuel, UPuZr):
            raise TypeError(
                f"fuel must be a UPuZr instance, got {type(fuel).__name__}"
            )
        if not isinstance(cladding, HT9):
            raise TypeError(
                f"cladding must be an HT9 instance, got {type(cladding).__name__}"
            )
        super().__init__(geometry)
        self._fuel     = fuel
        self._cladding = cladding

        # MNA-specific pending BCs
        self._gpres0       = None
        self._enable_zr    = None
        self._enable_waste = None
        self._enable_grsis = None
        self._sodium_mode  = None
        self._cond_model   = None
        self._Tpl_hist     = None

    def _make_cpp_solver(self):
        geom = self._geometry._make_cpp(_cpp)
        fuel = self._fuel._make_cpp(_cpp)
        clad = self._cladding._make_cpp(_cpp)
        return _cpp.FredMNaSolver(geom, fuel, clad)

    def _apply_bc_to_solver(self):
        super()._apply_bc_to_solver()
        if self._gpres0 is not None:
            self._solver.set_initial_gas_pressure(self._gpres0)
        if self._enable_zr is not None:
            self._solver.set_enable_zr_redistribution(self._enable_zr)
        if self._enable_waste is not None:
            self._solver.set_enable_clad_wastage(self._enable_waste)
        if self._enable_grsis is not None:
            self._solver.set_enable_grsis(self._enable_grsis)
        if self._sodium_mode is not None:
            mode = (_cpp.SodiumMode.TDependent if self._sodium_mode == 'TDependent'
                    else _cpp.SodiumMode.Constant)
            self._solver.set_sodium_mode(mode)
        if self._cond_model is not None:
            model_map = {
                'DetailedNaSodium': _cpp.ConductivityModel.DetailedNaSodium,
                'EmpiricalBurnup':  _cpp.ConductivityModel.EmpiricalBurnup,
                'EsfrSimple':       _cpp.ConductivityModel.EsfrSimple,
            }
            self._solver.set_conductivity_model(model_map[self._cond_model])
        if self._Tpl_hist is not None:
            self._solver.set_plenum_temperature_history(*self._Tpl_hist)

    # MNA-specific setters

    def set_initial_gas_pressure(self, p_MPa):
        """Initial sodium bond / fill-gas pressure [MPa]."""
        if p_MPa <= 0:
            raise ValueError(f"p_MPa must be positive, got {p_MPa!r}")
        self._gpres0 = p_MPa

    def set_plenum_temperature_history(self, times, T_K):
        from fred_solver_base import _check_timetable
        times, vals = _check_timetable("plenum temperature", times, T_K)
        for i, v in enumerate(vals):
            if v <= 0:
                raise ValueError(f"T_K[{i}] must be positive [K], got {v!r}")
        self._Tpl_hist = (times, vals)

    def set_enable_zr_redistribution(self, enable):
        """Enable/disable Zr redistribution model (default: enabled)."""
        self._enable_zr = bool(enable)

    def set_enable_clad_wastage(self, enable):
        """Enable/disable cladding wastage model (default: enabled)."""
        self._enable_waste = bool(enable)

    def set_enable_grsis(self, enable):
        """Enable/disable GRSIS bubble swelling model (default: enabled)."""
        self._enable_grsis = bool(enable)

    def set_sodium_mode(self, mode):
        """
        Sodium gap conductance mode.
        mode : 'TDependent' (default) or 'Constant'
        """
        if mode not in _SODIUM_MODES:
            raise ValueError(
                f"sodium_mode must be one of {sorted(_SODIUM_MODES)}, got {mode!r}"
            )
        self._sodium_mode = mode

    def set_conductivity_model(self, model):
        """
        Thermal conductivity correction model for U-Pu-Zr.
        model : 'DetailedNaSodium' (default), 'EmpiricalBurnup', 'EsfrSimple'
        """
        if model not in _CONDUCTIVITY_MODELS:
            raise ValueError(
                f"conductivity_model must be one of "
                f"{sorted(_CONDUCTIVITY_MODELS)}, got {model!r}"
            )
        self._cond_model = model

    # run

    def run(self, tend, dtout, output_file=None, all_steps=False, threads=1):
        """
        Run the simulation.

        Parameters
        ----------
        tend        : float  end time [s]
        dtout       : float  output interval [s]
        output_file : str    path to HDF5 output file; written incrementally
        all_steps   : bool   record every internal IDA step (debugging)
        threads     : int    number of OpenMP threads for the per-axial-layer
                              Newton solve (default 1 = serial).
        """
        if tend <= 0:
            raise ValueError(f"tend must be positive, got {tend!r}")
        if dtout <= 0:
            raise ValueError(f"dtout must be positive, got {dtout!r}")
        if dtout > tend:
            raise ValueError(f"dtout ({dtout}) must be <= tend ({tend})")
        if threads < 1:
            raise ValueError(f"threads must be >= 1, got {threads!r}")
        if self._channel is None:
            raise RuntimeError(
                "set_coolant_channel() must be called before run(). "
                "FRED-M-Na requires dhyd, xarea, flowr and an inlet temperature history."
            )

        self._build_solver()

        if output_file is not None:
            self._solver.set_output_file(str(output_file))
        self._solver.run(tend, dtout, bool(all_steps), int(threads))
        self._ran = True

    # Result accessors

    def time_points(self):
        self._require_run()
        return np.array(self._solver.times())

    def temperatures(self):
        """Temperature array [n_steps, nz, nf+nc] in K."""
        self._require_run()
        v = np.array(self._solver.temperatures())
        ns = len(self._solver.times())
        return v.reshape(ns, self._nz, self._nf + self._nc)

    def peak_fuel_temperature(self):
        self._require_run()
        t = self.temperatures()
        return t[:, :, 0].max(axis=1)

    def gas_pressure(self):
        self._require_run()
        return np.array(self._solver.gas_pressure())

    def fg_generated(self):
        self._require_run()
        return np.array(self._solver.fg_generated())

    def fg_released(self):
        self._require_run()
        return np.array(self._solver.fg_released())

    def gap_width(self):
        self._require_run()
        return np.array(self._solver.gap_width())

    def burnup(self):
        self._require_run()
        return np.array(self._solver.burnup())

    def clad_wastage(self):
        self._require_run()
        return np.array(self._solver.clad_wastage())

    def swelling_total(self):
        self._require_run()
        return np.array(self._solver.gris_swelling_total())

    def swelling_open(self):
        self._require_run()
        return np.array(self._solver.gris_swelling_open())

    def burst_margin(self):
        self._require_run()
        return np.array(self._solver.burst_margin())

    def melt_margin(self):
        self._require_run()
        return np.array(self._solver.melt_margin())


# ---------------------------------------------------------------------------
# read_results — read an HDF5 file written by FredMNaSolver.run()
# ---------------------------------------------------------------------------

def read_results(filename='results.h5'):
    """
    Read an HDF5 results file written by a FredMNaSolver run.

    Returns a dict with numpy arrays. Keys:
      time, T, peak_T_fuel, gap_width              — always present
      gpres, fggen, fgrel, bup                     — burnup/FGR quantities
      xwast, swtot, swopen, burst, melt            — MNA-specific
      restart_y, restart_yp                        — IDA state vectors
      nz, nf, nc, n_steps, app, all_steps          — metadata
    """
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py is required: conda install h5py")

    result = {}
    with h5py.File(filename, 'r') as f:
        meta = f['metadata']
        result['app']       = meta.attrs.get('app', 'fred-m-na')
        result['nz']        = int(meta.attrs['nz'])
        result['nf']        = int(meta.attrs['nf'])
        result['nc']        = int(meta.attrs['nc'])
        result['all_steps'] = bool(meta.attrs.get('all_steps', 0))

        result['time']    = f['time'][:]
        result['n_steps'] = len(result['time'])
        result['T']           = f['thermal/T'][:]
        result['peak_T_fuel'] = f['thermal/peak_T_fuel'][:]
        result['gap_width']   = f['thermal/gap_width'][:]

        if 'burnup' in f:
            bup = f['burnup']
            for key in ('gpres', 'fggen', 'fgrel', 'bup',
                        'xwast', 'swtot', 'swopen', 'burst', 'melt'):
                if key in bup:
                    result[key] = bup[key][:]

        if 'restart' in f:
            result['restart_y']  = f['restart/y'][:]
            result['restart_yp'] = f['restart/yp'][:]

    return result
