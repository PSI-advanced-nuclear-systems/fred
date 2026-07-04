"""
FRED-ROD — user-facing Python module.

Solves thermal conduction and stress-strain in a fuel rod with user-specified
material correlations. No irradiation physics is modelled.

Typical usage
-------------
    import sys; sys.path.insert(0, '/path/to/fred_platform/build')
    import fred_rod as fred

    g = fred.FuelRodGeometry()
    g.nf = 5; g.nc = 3; g.nz = 1
    g.rfi0 = 0.0;    g.rfo0 = 3.0e-3
    g.rci0 = 3.05e-3; g.rco0 = 3.4e-3
    g.dz0 = [0.10];  g.vgp = 1e-6
    g.build()

    fuel   = fred.UO2()
    clad   = fred.AIM1()
    gap    = fred.HeliumGap()
    solver = fred.FredRodSolver(g, fuel, clad, gap)

    solver.set_initial_temperature(900.0)
    solver.set_coolant_temperature_history([0, 1e5], [773.0, 773.0])
    solver.set_coolant_pressure_history([0, 1e5], [0.1, 0.1])
    solver.set_power_density_history([0, 1e5], [3e8, 3e8])
    solver.run(tend=1e5, dtout=1e4, output_file='results.h5')

    r = fred.read_results('results.h5')
"""

import numpy as np

import _fred_rod as _cpp   # compiled C++ extension

from fred_geometry import FuelRodGeometry
from fred_materials import FuelPellet, Cladding, GapFill, _validate_density
from fred_solver_base import FredSolverBasePrescribedCoolant


# ---------------------------------------------------------------------------
# Material wrappers
# ---------------------------------------------------------------------------

class UO2(FuelPellet):
    """UO2 fuel pellet (Fink 2000 correlations)."""
    def __init__(self, reference_density=None):
        _validate_density("UO2", reference_density)
        self._ref_dens = reference_density
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = (cpp_module.UO2(self._ref_dens)
                             if self._ref_dens is not None else cpp_module.UO2())
        return self._cpp_obj


class MOX(FuelPellet):
    """MOX fuel pellet. pu_content: Pu weight fraction (e.g. 0.20 for 20 wt%)."""
    def __init__(self, pu_content, reference_density=None):
        if not (0 < pu_content < 1):
            raise ValueError(f"pu_content must be in (0, 1), got {pu_content!r}")
        _validate_density("MOX", reference_density)
        self._pu = pu_content
        self._ref_dens = reference_density
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = (cpp_module.MOX(self._pu, self._ref_dens)
                             if self._ref_dens is not None else cpp_module.MOX(self._pu))
        return self._cpp_obj


class DummyFuelPellet(FuelPellet):
    """Placeholder fuel pellet for geometry-only or test runs."""
    def __init__(self):
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = cpp_module.DummyFuelPellet()
        return self._cpp_obj


class AIM1(Cladding):
    """AIM-1 austenitic steel cladding."""
    def __init__(self, reference_density=None):
        _validate_density("AIM1", reference_density)
        self._ref_dens = reference_density
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = (cpp_module.AIM1(self._ref_dens)
                             if self._ref_dens is not None else cpp_module.AIM1())
        return self._cpp_obj


class T91(Cladding):
    """T91 ferritic-martensitic steel cladding."""
    def __init__(self, reference_density=None):
        _validate_density("T91", reference_density)
        self._ref_dens = reference_density
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = (cpp_module.T91(self._ref_dens)
                             if self._ref_dens is not None else cpp_module.T91())
        return self._cpp_obj


class DummyCladding(Cladding):
    """Placeholder cladding for geometry-only or test runs."""
    def __init__(self):
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = cpp_module.DummyCladding()
        return self._cpp_obj


class HeliumGap(GapFill):
    """Pure helium gap fill gas."""
    def __init__(self):
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = cpp_module.He()
        return self._cpp_obj


class HeliumKryptonXenonGap(GapFill):
    """He-Kr-Xe mixed fill gas."""
    def __init__(self):
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = cpp_module.HeKrXe()
        return self._cpp_obj


class DummyGapMaterial(GapFill):
    """Placeholder gap material."""
    def __init__(self):
        self._cpp_obj = None

    def _make_cpp(self, cpp_module):
        if self._cpp_obj is None:
            self._cpp_obj = cpp_module.DummyGapMaterial()
        return self._cpp_obj


# ---------------------------------------------------------------------------
# FredRodSolver
# ---------------------------------------------------------------------------

class FredRodSolver(FredSolverBasePrescribedCoolant):
    """
    FRED-ROD solver: thermal conduction + stress-strain.
    No irradiation physics — user specifies material correlations directly.
    """

    def __init__(self, geometry, fuel, cladding, gap):
        if not isinstance(fuel, FuelPellet):
            raise TypeError(
                f"fuel must be a FuelPellet instance (e.g. UO2(), MOX()), "
                f"got {type(fuel).__name__}"
            )
        if not isinstance(cladding, Cladding):
            raise TypeError(
                f"cladding must be a Cladding instance (e.g. AIM1(), T91()), "
                f"got {type(cladding).__name__}"
            )
        if not isinstance(gap, GapFill):
            raise TypeError(
                f"gap must be a GapFill instance (e.g. HeliumGap()), "
                f"got {type(gap).__name__}"
            )
        super().__init__(geometry)
        self._fuel    = fuel
        self._cladding = cladding
        self._gap     = gap

        # ROD-specific pending BCs
        self._gpres0  = None
        self._heat_on = None
        self._mech_on = None

    def _make_cpp_solver(self):
        geom = self._geometry._make_cpp(_cpp)
        fuel = self._fuel._make_cpp(_cpp)
        clad = self._cladding._make_cpp(_cpp)
        gap  = self._gap._make_cpp(_cpp)
        return _cpp.FredRodSolver(geom, fuel, clad, gap)

    def _apply_bc_to_solver(self):
        super()._apply_bc_to_solver()
        if self._gpres0 is not None:
            self._solver.set_initial_gas_pressure(self._gpres0)
        if self._heat_on is not None:
            self._solver.set_enable_heat_conduction(self._heat_on)
        if self._mech_on is not None:
            self._solver.set_enable_stress_strain(self._mech_on)

    # ROD-specific setters

    def set_initial_gas_pressure(self, gpres0_MPa):
        """Initial fill-gas pressure [MPa]. Must be called before run()."""
        if gpres0_MPa <= 0:
            raise ValueError(f"gpres0_MPa must be positive, got {gpres0_MPa!r}")
        self._gpres0 = gpres0_MPa

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
                              residual loop (default 1 = serial). Silently
                              clamped to 1 if a Python-defined material
                              subclass is in use (not thread-safe).
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
        stride = self._nz * (self._nf + self._nc)
        return v.reshape(ns, self._nz, self._nf + self._nc) if v.size == ns * stride else v

    def peak_fuel_temperature(self):
        self._require_run()
        return np.array(self._solver.peak_fuel_temperature())

    def gap_width(self):
        self._require_run()
        return np.array(self._solver.gap_width())

    def clad_outer_hoop_stress(self):
        self._require_run()
        return np.array(self._solver.clad_outer_hoop_stress())

    def clad_outer_radius(self):
        self._require_run()
        return np.array(self._solver.clad_outer_radius())

    def fuel_outer_radius(self):
        self._require_run()
        return np.array(self._solver.fuel_outer_radius())

    def clad_inner_radius(self):
        self._require_run()
        return np.array(self._solver.clad_inner_radius())

    def contact_pressure(self):
        self._require_run()
        return np.array(self._solver.contact_pressure())


# ---------------------------------------------------------------------------
# read_results — read an HDF5 file written by FredRodSolver.run()
# ---------------------------------------------------------------------------

def read_results(filename='results.h5'):
    """
    Read an HDF5 results file written by a FredRodSolver run.

    Returns a dict with numpy arrays. Keys:
      time, T, peak_T_fuel, gap_width       — always present
      pfc, rfo, rci, rco, sigh_outer        — mechanical (stress-strain on)
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
        result['app']       = meta.attrs.get('app', 'fred-rod')
        result['nz']        = int(meta.attrs['nz'])
        result['nf']        = int(meta.attrs['nf'])
        result['nc']        = int(meta.attrs['nc'])
        result['all_steps'] = bool(meta.attrs.get('all_steps', 0))

        result['time']    = f['time'][:]
        result['n_steps'] = len(result['time'])
        result['T']           = f['thermal/T'][:]
        result['peak_T_fuel'] = f['thermal/peak_T_fuel'][:]
        result['gap_width']   = f['thermal/gap_width'][:]

        if 'mechanical' in f:
            for key in ('pfc', 'rfo', 'rci', 'rco', 'sigh_outer'):
                result[key] = f[f'mechanical/{key}'][:]

        if 'restart' in f:
            result['restart_y']  = f['restart/y'][:]
            result['restart_yp'] = f['restart/yp'][:]

    return result
