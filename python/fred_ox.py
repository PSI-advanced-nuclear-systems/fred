"""
FRED-OX — user-facing Python module.

Extends FRED-ROD with MOX fuel physics: base irradiation effects on thermal
conductivity, empirical fission gas production/release, and fuel swelling.

Typical usage
-------------
    import sys; sys.path.insert(0, '/path/to/fred_platform/build')
    import fred_ox as fred

    g = fred.FuelRodGeometry()
    g.nf = 5; g.nc = 3; g.nz = 1
    g.rfi0 = 0.0;    g.rfo0 = 3.0e-3
    g.rci0 = 3.05e-3; g.rco0 = 3.4e-3
    g.dz0 = [0.10];  g.vgp = 1e-6
    g.build()

    fuel   = fred.FredOxMOX(pu_content=0.30)
    clad   = fred.FredOxAIM1()
    gap    = fred.FredOxGapMaterial()
    solver = fred.FredOxSolver(g, fuel, clad, gap)

    solver.set_initial_temperature(900.0)
    solver.set_coolant_temperature_history([0, 1e7], [723.0, 723.0])
    solver.set_coolant_pressure_history([0, 1e7], [0.1, 0.1])
    solver.set_power_density_history([0, 1e7], [3e8, 3e8])
    solver.set_initial_gas_pressure(0.1)
    solver.run(tend=1e7, dtout=1e6, output_file='results.h5')

    r = fred.read_results('results.h5')
"""

import numpy as np

import _fred_ox as _cpp   # compiled C++ extension

from fred_geometry import FuelRodGeometry
from fred_materials import FuelPellet, Cladding, GapFill, _validate_density
from fred_solver_base import FredSolverBasePrescribedCoolant


# ---------------------------------------------------------------------------
# OX-specific material wrappers
# ---------------------------------------------------------------------------

class FredOxMOX(FuelPellet):
    """
    MOX fuel pellet with burnup-dependent thermal conductivity (FRED-OX).

    Parameters
    ----------
    pu_content              : Pu weight fraction [-], e.g. 0.30 for 30 wt%
    reference_density_frac  : as-fabricated density as fraction of theoretical
                              density [-] (default 1.0)
    stoichiometry           : O/M ratio (default 2.0, range 1.9–2.0)
    """
    def __init__(self, pu_content, reference_density_frac=1.0, stoichiometry=2.0):
        if not (0 < pu_content < 1):
            raise ValueError(f"pu_content must be in (0, 1), got {pu_content!r}")
        if not (0 < reference_density_frac <= 1):
            raise ValueError(
                f"reference_density_frac must be in (0, 1], got {reference_density_frac!r}"
            )
        if not (1.9 <= stoichiometry <= 2.0):
            raise ValueError(
                f"stoichiometry (O/M) must be in [1.9, 2.0], got {stoichiometry!r}"
            )
        self._pu  = pu_content
        self._rof = reference_density_frac
        self._sto = stoichiometry
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            # The C++ constructor's second argument is an *absolute* density
            # [kg/m3], not the [0, 1] theoretical-density fraction this
            # wrapper takes — convert using the same theoretical-density
            # correlation the C++ side uses internally (via a throwaway
            # instance, since theoretical density only depends on pu_content).
            tden = cpp_module.FredOxMOX(self._pu, 1.0, self._sto).theoretical_density()
            abs_density = self._rof * tden
            self._cpp_obj = cpp_module.FredOxMOX(self._pu, abs_density, self._sto)
        return self._cpp_obj


class FredOxAIM1(Cladding):
    """AIM-1 cladding with burnup-dependent properties (FRED-OX)."""
    def __init__(self, reference_density=None):
        _validate_density("FredOxAIM1", reference_density)
        self._ref_dens = reference_density
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = (cpp_module.FredOxAIM1(self._ref_dens)
                             if self._ref_dens is not None else cpp_module.FredOxAIM1())
        return self._cpp_obj


class FredOxT91(Cladding):
    """T91 cladding with burnup-dependent properties (FRED-OX)."""
    def __init__(self, reference_density=None):
        _validate_density("FredOxT91", reference_density)
        self._ref_dens = reference_density
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = (cpp_module.FredOxT91(self._ref_dens)
                             if self._ref_dens is not None else cpp_module.FredOxT91())
        return self._cpp_obj


class FredOxGapMaterial(GapFill):
    """Gap material for FRED-OX (mixed fission gas conductivity model)."""
    def __init__(self):
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = cpp_module.FredOxGapMaterial()
        return self._cpp_obj


# ---------------------------------------------------------------------------
# FredOxSolver
# ---------------------------------------------------------------------------

class FredOxSolver(FredSolverBasePrescribedCoolant):
    """
    FRED-OX solver: thermal conduction + stress-strain + MOX irradiation physics.
    Includes fission gas production/release and fuel swelling (MATPRO).
    """

    def __init__(self, geometry, mox, cladding, gap):
        if not isinstance(mox, FredOxMOX):
            raise TypeError(
                f"mox must be a FredOxMOX instance, got {type(mox).__name__}"
            )
        if not isinstance(cladding, Cladding):
            raise TypeError(
                f"cladding must be a Cladding instance (e.g. FredOxAIM1()), "
                f"got {type(cladding).__name__}"
            )
        if not isinstance(gap, GapFill):
            raise TypeError(
                f"gap must be a GapFill instance (e.g. FredOxGapMaterial()), "
                f"got {type(gap).__name__}"
            )
        super().__init__(geometry)
        self._mox      = mox
        self._cladding = cladding
        self._gap      = gap

        # OX-specific pending BCs
        self._gpres0      = None
        self._fswelmlt    = None
        self._heat_on     = None
        self._mech_on     = None

    def _make_cpp_solver(self):
        geom = self._geometry._make_cpp(_cpp)
        mox  = self._mox._make_cpp(_cpp)
        clad = self._cladding._make_cpp(_cpp)
        gap  = self._gap._make_cpp(_cpp)
        return _cpp.FredOxSolver(geom, mox, clad, gap)

    def _apply_bc_to_solver(self):
        super()._apply_bc_to_solver()
        if self._gpres0 is not None:
            self._solver.set_initial_gas_pressure(self._gpres0)
        if self._fswelmlt is not None:
            self._solver.set_swelling_multiplier(self._fswelmlt)
        if self._heat_on is not None:
            self._solver.set_enable_heat_conduction(self._heat_on)
        if self._mech_on is not None:
            self._solver.set_enable_stress_strain(self._mech_on)

    # OX-specific setters

    def set_initial_gas_pressure(self, gpres0_MPa):
        """Initial fill-gas pressure [MPa]. Must be called before run()."""
        if gpres0_MPa <= 0:
            raise ValueError(f"gpres0_MPa must be positive, got {gpres0_MPa!r}")
        self._gpres0 = gpres0_MPa

    def set_swelling_multiplier(self, fswelmlt):
        """Fuel swelling multiplier (1.0 = full MATPRO, 0.8 = legacy default)."""
        if not (0 < fswelmlt <= 2):
            raise ValueError(f"fswelmlt must be in (0, 2], got {fswelmlt!r}")
        self._fswelmlt = fswelmlt

    def set_enable_heat_conduction(self, enable):
        self._heat_on = bool(enable)

    def set_enable_stress_strain(self, enable):
        self._mech_on = bool(enable)

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
                              residual loop (default 1 = serial).
        """
        if tend <= 0:
            raise ValueError(f"tend must be positive, got {tend!r}")
        if dtout <= 0:
            raise ValueError(f"dtout must be positive, got {dtout!r}")
        if dtout > tend:
            raise ValueError(f"dtout ({dtout}) must be <= tend ({tend})")
        if threads < 1:
            raise ValueError(f"threads must be >= 1, got {threads!r}")

        self._build_solver()

        if output_file is not None:
            self._solver.set_output_file(str(output_file))
        self._solver.run(tend, dtout, bool(all_steps), int(threads))
        self._ran = True

    # Result accessors

    def time_points(self):
        self._require_run()
        return np.array(self._solver.time_points())

    def temperatures(self):
        """Temperature array [n_steps, nz, nf+nc] in K."""
        self._require_run()
        v = np.array(self._solver.temperatures())
        ns = len(self._solver.time_points())
        return v.reshape(ns, self._nz, self._nf + self._nc)

    def peak_fuel_temperature(self):
        self._require_run()
        return np.array(self._solver.peak_fuel_temperature())

    def gap_width(self):
        self._require_run()
        return np.array(self._solver.gap_width())

    def gas_pressure(self):
        self._require_run()
        return np.array(self._solver.gas_pressure())

    def fg_generated(self):
        self._require_run()
        return np.array(self._solver.fg_generated())

    def fg_released(self):
        self._require_run()
        return np.array(self._solver.fg_released())

    def burnup(self):
        self._require_run()
        return np.array(self._solver.burnup())


# ---------------------------------------------------------------------------
# read_results — read an HDF5 file written by FredOxSolver.run()
# ---------------------------------------------------------------------------

def read_results(filename='results.h5'):
    """
    Read an HDF5 results file written by a FredOxSolver run.

    Returns a dict with numpy arrays. Keys:
      time, T, peak_T_fuel, gap_width       — always present
      gpres, fggen, fgrel, bup              — burnup/FGR quantities
      restart_y, restart_yp                 — IDA state vectors
      nz, nf, nc, n_steps, app, all_steps   — metadata
    """
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py is required: conda install h5py")

    result = {}
    with h5py.File(filename, 'r') as f:
        meta = f['metadata']
        result['app']       = meta.attrs.get('app', 'fred-ox')
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
            for key in ('gpres', 'fggen', 'fgrel', 'bup'):
                result[key] = f[f'burnup/{key}'][:]

        if 'restart' in f:
            result['restart_y']  = f['restart/y'][:]
            result['restart_yp'] = f['restart/yp'][:]

    return result
