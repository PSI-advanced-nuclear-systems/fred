"""
Python solver base class hierarchy for the FRED platform.

FredSolverBase
├── FredSolverBasePrescribedCoolant  (FRED-ROD, FRED-OX)
└── FredSolverBaseSubchannel         (FRED-M-Na)
"""

import numpy as np


# ---------------------------------------------------------------------------
# Shared validation helpers
# ---------------------------------------------------------------------------

def _check_timetable(name, times, values):
    times = list(times)
    values = list(values)
    if len(times) < 2:
        raise ValueError(f"{name}: at least 2 time points required, got {len(times)}")
    for i in range(1, len(times)):
        if times[i] <= times[i - 1]:
            raise ValueError(
                f"{name}: times must be strictly increasing; "
                f"got {times[i-1]} then {times[i]} at index {i}"
            )
    if len(times) != len(values):
        raise ValueError(
            f"{name}: times has {len(times)} entries but values has {len(values)}"
        )
    return times, values


def _check_per_layer_timetable(name, times, per_layer, nz):
    times = list(times)
    if len(times) < 2:
        raise ValueError(f"{name}: at least 2 time points required, got {len(times)}")
    for i in range(1, len(times)):
        if times[i] <= times[i - 1]:
            raise ValueError(
                f"{name}: times must be strictly increasing; "
                f"got {times[i-1]} then {times[i]} at index {i}"
            )
    if len(per_layer) != nz:
        raise ValueError(
            f"{name}: expected {nz} value lists (one per axial layer), "
            f"got {len(per_layer)}"
        )
    for j, vals in enumerate(per_layer):
        if len(vals) != len(times):
            raise ValueError(
                f"{name}: layer {j} has {len(vals)} values but "
                f"{len(times)} time points were supplied"
            )
    return times, [list(v) for v in per_layer]


# ---------------------------------------------------------------------------
# FredSolverBase — common to all apps
# ---------------------------------------------------------------------------

class FredSolverBase:
    """
    Common base for all FRED solver wrappers.
    Subclasses must implement run() and _make_cpp_solver().
    """

    def __init__(self, geometry):
        from fred_geometry import FuelRodGeometry
        if not isinstance(geometry, FuelRodGeometry):
            raise TypeError(
                f"geometry must be a FuelRodGeometry instance, got {type(geometry).__name__}"
            )
        if not geometry._built:
            geometry.build()
        self._geometry = geometry
        self._nz = geometry.nz
        self._nf = geometry.nf
        self._nc = geometry.nc
        self._solver = None
        self._ran = False
        self._output_file = None

        # Pending BC state (applied to C++ solver after construction)
        self._T0        = None
        self._qqv_hist  = None      # (times, values) uniform
        self._qqv_layer = None      # (times, per_layer)
        self._rtol      = None
        self._atol      = None

    # -- common setters ------------------------------------------------------

    def set_initial_temperature(self, T0_K):
        if T0_K <= 0:
            raise ValueError(f"T0_K must be positive [K], got {T0_K!r}")
        self._T0 = T0_K
        if self._solver is not None:
            self._solver.set_initial_temperature(T0_K)

    def set_power_density_history(self, times, qqv_W_m3):
        times, vals = _check_timetable("power density", times, qqv_W_m3)
        for i, v in enumerate(vals):
            if v < 0:
                raise ValueError(f"qqv_W_m3[{i}] must be >= 0, got {v!r}")
        self._qqv_hist = (times, vals)
        self._qqv_layer = None
        if self._solver is not None:
            self._solver.set_power_density_history(times, vals)

    def set_power_density_history_per_layer(self, times, qqv_per_layer):
        times, layers = _check_per_layer_timetable(
            "power density per layer", times, qqv_per_layer, self._nz)
        self._qqv_layer = (times, layers)
        self._qqv_hist = None
        if self._solver is not None:
            self._solver.set_power_density_history_per_layer(times, layers)

    def set_tolerances(self, rtol, atol):
        if rtol <= 0 or atol <= 0:
            raise ValueError(f"rtol and atol must be positive, got rtol={rtol!r}, atol={atol!r}")
        self._rtol = rtol
        self._atol = atol
        if self._solver is not None:
            self._solver.set_tolerances(rtol, atol)

    def set_output_file(self, path):
        self._output_file = str(path)

    # -- lifecycle -----------------------------------------------------------

    def _build_solver(self):
        """Construct the C++ solver and apply all pending BCs. Called by run()."""
        self._solver = self._make_cpp_solver()
        if self._T0 is not None:
            self._solver.set_initial_temperature(self._T0)
        if self._qqv_hist is not None:
            self._solver.set_power_density_history(*self._qqv_hist)
        if self._qqv_layer is not None:
            self._solver.set_power_density_history_per_layer(*self._qqv_layer)
        if self._rtol is not None:
            self._solver.set_tolerances(self._rtol, self._atol)
        self._apply_bc_to_solver()

    def _apply_bc_to_solver(self):
        """Apply subclass-specific BCs to self._solver. Override in subclasses."""

    def _make_cpp_solver(self):
        """Create and return the C++ solver object. Must be implemented by subclasses."""
        raise NotImplementedError

    def _require_run(self):
        if not self._ran:
            raise RuntimeError("call run() before accessing results")

    def run(self, tend, dtout, output_file=None, all_steps=False, threads=1):
        raise NotImplementedError


# ---------------------------------------------------------------------------
# FredSolverBasePrescribedCoolant — coolant T and P given as time tables
# ---------------------------------------------------------------------------

class FredSolverBasePrescribedCoolant(FredSolverBase):
    """
    Base for solvers where the user prescribes coolant temperature and pressure
    directly as time-varying boundary conditions (FRED-ROD, FRED-OX).
    """

    def __init__(self, geometry):
        super().__init__(geometry)
        self._Tc_hist   = None   # (times, T_K) uniform
        self._Tc_layer  = None   # (times, per_layer)
        self._pc_hist   = None   # (times, pcool_MPa)
        self._Tpl_hist  = None   # (times, T_K) plenum

    def set_coolant_temperature_history(self, times, T_K):
        times, vals = _check_timetable("coolant temperature", times, T_K)
        for i, v in enumerate(vals):
            if v <= 0:
                raise ValueError(f"T_K[{i}] must be positive [K], got {v!r}")
        self._Tc_hist  = (times, vals)
        self._Tc_layer = None
        if self._solver is not None:
            self._solver.set_coolant_temperature(times, vals)

    def set_coolant_temperature_history_per_layer(self, times, T_per_layer):
        times, layers = _check_per_layer_timetable(
            "coolant temperature per layer", times, T_per_layer, self._nz)
        for j, vals in enumerate(layers):
            for i, v in enumerate(vals):
                if v <= 0:
                    raise ValueError(
                        f"coolant temperature per layer: layer {j}, entry {i} "
                        f"must be positive [K], got {v!r}"
                    )
        self._Tc_layer = (times, layers)
        self._Tc_hist  = None
        if self._solver is not None:
            self._solver.set_coolant_temperature_per_layer(times, layers)

    def set_coolant_pressure_history(self, times, pcool_MPa):
        times, vals = _check_timetable("coolant pressure", times, pcool_MPa)
        for i, v in enumerate(vals):
            if v <= 0:
                raise ValueError(f"pcool_MPa[{i}] must be positive, got {v!r}")
        self._pc_hist = (times, vals)
        if self._solver is not None:
            self._solver.set_coolant_pressure_history(times, vals)

    def set_plenum_temperature_history(self, times, T_K):
        times, vals = _check_timetable("plenum temperature", times, T_K)
        for i, v in enumerate(vals):
            if v <= 0:
                raise ValueError(f"T_K[{i}] must be positive [K], got {v!r}")
        self._Tpl_hist = (times, vals)
        if self._solver is not None:
            self._solver.set_plenum_temperature_history(times, vals)

    def _apply_bc_to_solver(self):
        if self._Tc_hist is not None:
            self._solver.set_coolant_temperature(*self._Tc_hist)
        if self._Tc_layer is not None:
            self._solver.set_coolant_temperature_per_layer(*self._Tc_layer)
        if self._pc_hist is not None:
            self._solver.set_coolant_pressure_history(*self._pc_hist)
        if self._Tpl_hist is not None:
            self._solver.set_plenum_temperature_history(*self._Tpl_hist)


# ---------------------------------------------------------------------------
# FredSolverBaseSubchannel — coolant computed from subchannel energy balance
# ---------------------------------------------------------------------------

class FredSolverBaseSubchannel(FredSolverBase):
    """
    Base for solvers that compute coolant conditions from a subchannel energy
    balance with a Robin BC at the cladding outer surface (FRED-M-Na).
    Prescribing coolant temperature directly is not supported.
    """

    def __init__(self, geometry):
        super().__init__(geometry)
        self._channel = None   # dict of set_coolant_channel args
        self._pc_hist = None   # plenum/system gas pressure history

    def set_coolant_channel(self, dhyd, xarea, flowr,
                            T_inlet_times, T_inlet_vals,
                            htc_correlation='Mikityuk'):
        _CORR = {'Mikityuk', 'Subbotin'}
        if dhyd <= 0:
            raise ValueError(f"dhyd must be positive [m], got {dhyd!r}")
        if xarea <= 0:
            raise ValueError(f"xarea must be positive [m²], got {xarea!r}")
        if flowr <= 0:
            raise ValueError(f"flowr must be positive [kg/s], got {flowr!r}")
        T_inlet_times, T_inlet_vals = _check_timetable(
            "inlet temperature", T_inlet_times, T_inlet_vals)
        for i, v in enumerate(T_inlet_vals):
            if v <= 0:
                raise ValueError(
                    f"T_inlet_vals[{i}] must be positive [K], got {v!r}"
                )
        if htc_correlation not in _CORR:
            raise ValueError(
                f"htc_correlation must be one of {sorted(_CORR)}, "
                f"got {htc_correlation!r}"
            )
        self._channel = dict(
            dhyd=dhyd, xarea=xarea, flowr=flowr,
            times=T_inlet_times, vals=T_inlet_vals,
            corr=htc_correlation,
        )
        if self._solver is not None:
            self._apply_channel()

    def set_coolant_pressure_history(self, times, pcool_MPa):
        times, vals = _check_timetable("coolant pressure", times, pcool_MPa)
        for i, v in enumerate(vals):
            if v <= 0:
                raise ValueError(f"pcool_MPa[{i}] must be positive, got {v!r}")
        self._pc_hist = (times, vals)
        if self._solver is not None:
            self._solver.set_coolant_pressure_history(*self._pc_hist)

    def set_coolant_temperature_history(self, *args, **kwargs):
        raise AttributeError(
            "This solver uses subchannel calculation for coolant temperature.\n"
            "Use set_coolant_channel(dhyd, xarea, flowr, T_inlet_times, T_inlet_vals) instead."
        )

    def _apply_channel(self):
        if self._channel is None:
            return
        c = self._channel
        self._solver.set_coolant_channel(
            c['dhyd'], c['xarea'], c['flowr'],
            c['times'], c['vals'], c['corr'],
        )

    def _apply_bc_to_solver(self):
        self._apply_channel()
        if self._pc_hist is not None:
            self._solver.set_coolant_pressure_history(*self._pc_hist)
