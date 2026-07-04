"""
fred_input — validated Python wrapper for the FRED 2.0 platform API.

Replace ``import fred`` with ``import fred_input as fred`` (or
``from fred_input import *``) for fully validated input.  All geometry,
material, and solver parameters are checked before being forwarded to the
C++ layer, so errors are reported in Python with clear, plain-English
messages.

Why validate in Python?
  - Errors surface immediately, before any C++ object is constructed.
  - Validation rules can be extended without recompiling.
  - The C++ layer can trust its inputs and use assert() only for internal
    invariants.

Usage (drop-in replacement for bare ``fred``):

    import sys
    sys.path.insert(0, '/path/to/fred_platform/build')
    import fred_input as fred          # validated API

    g = fred.FuelRodGeometry()
    g.nf = 3; g.nc = 2; g.nz = 1
    g.rfi0 = 0.0; g.rfo0 = 4.2e-3
    g.rci0 = 4.215e-3; g.rco0 = 4.9e-3
    g.dz0 = [0.10]; g.vgp = 1e-6
    g.ruff = 5e-6; g.rufc = 5e-6
    g.build()

    fuel   = fred.DummyFuelPellet()
    clad   = fred.DummyCladding()
    gap    = fred.DummyGapMaterial()
    solver = fred.FredRodSolver(g, fuel, clad, gap)
    solver.set_power_density_history([0, 1e6], [1e8, 1e8])
    solver.set_coolant_temperature([0, 1e6], [700.0, 700.0])
    solver.set_coolant_pressure(5.0)
    solver.set_initial_temperature(700.0)
    solver.set_initial_gas_pressure(0.1)
    solver.run(tend=10.0, dtout=1.0)
"""

import fred as _fred  # the compiled C++ extension

# Re-export abstract base classes unchanged (users need them for custom materials).
FuelPelletMaterial = _fred.FuelPelletMaterial
CladdingMaterial   = _fred.CladdingMaterial
GapMaterial        = _fred.GapMaterial

# ---------------------------------------------------------------------------
# Internal validation helpers
# ---------------------------------------------------------------------------

def _pos(name, v):
    """Require v > 0."""
    if v <= 0:
        raise ValueError(f"{name} must be positive, got {v!r}")


def _nonneg(name, v):
    """Require v >= 0."""
    if v < 0:
        raise ValueError(f"{name} must be non-negative, got {v!r}")


def _int_ge(name, v, minimum):
    """Require integer v >= minimum."""
    if not isinstance(v, int):
        raise TypeError(f"{name} must be an integer, got {type(v).__name__!r}")
    if v < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {v}")


def _range_incl(name, v, lo, hi):
    """Require lo <= v <= hi."""
    if not (lo <= v <= hi):
        raise ValueError(f"{name} must be in [{lo}, {hi}], got {v!r}")


def _monotone_times(label, times):
    """Require a strictly increasing sequence with at least 2 entries."""
    if len(times) < 2:
        raise ValueError(
            f"{label}: at least 2 time points required, got {len(times)}"
        )
    for i in range(1, len(times)):
        if times[i] <= times[i - 1]:
            raise ValueError(
                f"{label}: time points must be strictly increasing; "
                f"got {times[i - 1]} then {times[i]} at index {i}"
            )


def _same_len(label_a, a, label_b, b):
    if len(a) != len(b):
        raise ValueError(
            f"{label_a} has {len(a)} entries but {label_b} has {len(b)}"
        )


def _all_nonneg_vals(label, vals):
    for i, v in enumerate(vals):
        if v < 0:
            raise ValueError(
                f"{label}[{i}] must be non-negative, got {v!r}"
            )


def _all_positive_K(label, vals):
    for i, v in enumerate(vals):
        if v <= 0:
            raise ValueError(
                f"{label}[{i}] must be a positive temperature (K), got {v!r}"
            )


def _all_positive_MPa(label, vals):
    for i, v in enumerate(vals):
        if v <= 0:
            raise ValueError(
                f"{label}[{i}] must be positive (MPa), got {v!r}"
            )


def _check_time_history(label, times, values, check_values_fn=None):
    """Common validation for a (times, values) boundary-condition history."""
    _monotone_times(label + " times", times)
    _same_len(label + " times", times, label + " values", values)
    if check_values_fn is not None:
        check_values_fn(label + " values", values)


def _check_per_layer_history(label, times, per_layer, nz):
    """Validate a per-layer BC history: per_layer[j] is the value list for layer j."""
    _monotone_times(label + " times", times)
    if len(per_layer) != nz:
        raise ValueError(
            f"{label}: expected {nz} layer lists (one per axial layer), "
            f"got {len(per_layer)}"
        )
    nt = len(times)
    for j, vals in enumerate(per_layer):
        if len(vals) != nt:
            raise ValueError(
                f"{label}: layer {j} has {len(vals)} values but "
                f"{nt} time points were supplied"
            )


# ---------------------------------------------------------------------------
# FuelRodGeometry
# ---------------------------------------------------------------------------

class FuelRodGeometry:
    """
    Validated Python wrapper for fred.FuelRodGeometry.

    Set all named attributes, then call build().  After build() the
    underlying C++ object is available as ._cpp; additional derived
    attributes (drf0, drc0, rad0, …) are forwarded to the C++ object.

    Required attributes
    -------------------
    nf   : int >= 2          number of fuel radial nodes
    nc   : int >= 2          number of cladding radial nodes
    nz   : int >= 1          number of axial layers
    rfi0 : float >= 0        inner fuel radius [m]  (0 for solid pellet)
    rfo0 : float > rfi0      outer fuel radius [m]
    rci0 : float > rfo0      inner cladding radius [m]  (gap must be positive)
    rco0 : float > rci0      outer cladding radius [m]
    dz0  : list of nz floats axial layer heights [m]  (all > 0)

    Optional attributes (have sensible defaults)
    --------------------------------------------
    vgp  : float >= 0        gas plenum volume [m3]  (default 0.0)
    ruff : float > 0         fuel outer surface roughness [m]  (default 1 µm)
    rufc : float > 0         cladding inner surface roughness [m]  (default 1 µm)
    """

    def __init__(self):
        # Settable parameters (None = not yet assigned)
        self.nf   = None
        self.nc   = None
        self.nz   = None
        self.rfi0 = None
        self.rfo0 = None
        self.rci0 = None
        self.rco0 = None
        self.dz0  = None
        # Defaults matching C++ struct defaults
        self.vgp  = 0.0
        self.ruff = 1.0e-6
        self.rufc = 1.0e-6
        self._cpp_obj = None

    # -- transparent attribute forwarding after build() -------------------
    def __getattr__(self, name):
        # Only called when normal attribute lookup fails (i.e. not in __dict__)
        if name.startswith('_'):
            raise AttributeError(name)
        cpp = self.__dict__.get('_cpp_obj')
        if cpp is not None:
            return getattr(cpp, name)
        raise AttributeError(
            f"FuelRodGeometry has no attribute {name!r}; "
            "did you forget to call build()?"
        )

    # -- validation + construction ----------------------------------------
    def build(self):
        """
        Validate all geometry parameters and construct the C++ object.

        Raises ValueError / TypeError with a plain-English message for any
        invalid combination of parameters.
        """
        # Node counts
        for attr in ("nf", "nc", "nz"):
            if getattr(self, attr) is None:
                raise ValueError(f"{attr} is not set")
        _int_ge("nf", self.nf, 2)
        _int_ge("nc", self.nc, 2)
        _int_ge("nz", self.nz, 1)

        # Radii
        for attr in ("rfi0", "rfo0", "rci0", "rco0"):
            if getattr(self, attr) is None:
                raise ValueError(f"{attr} is not set")
        _nonneg("rfi0", self.rfi0)
        _pos("rfo0", self.rfo0)
        if self.rfo0 <= self.rfi0:
            raise ValueError(
                f"rfo0 ({self.rfo0:.6g} m) must be greater than "
                f"rfi0 ({self.rfi0:.6g} m)"
            )
        if self.rci0 <= self.rfo0:
            raise ValueError(
                f"rci0 ({self.rci0:.6g} m) must be greater than "
                f"rfo0 ({self.rfo0:.6g} m) — as-fabricated gap must be positive"
            )
        if self.rco0 <= self.rci0:
            raise ValueError(
                f"rco0 ({self.rco0:.6g} m) must be greater than "
                f"rci0 ({self.rci0:.6g} m)"
            )

        # Axial layers
        if self.dz0 is None:
            raise ValueError("dz0 is not set")
        dz0 = list(self.dz0)
        if len(dz0) != self.nz:
            raise ValueError(
                f"dz0 must have nz={self.nz} entries, got {len(dz0)}"
            )
        for i, dz in enumerate(dz0):
            if dz <= 0:
                raise ValueError(
                    f"dz0[{i}] must be positive [m], got {dz!r}"
                )

        # Auxiliary
        _nonneg("vgp",  self.vgp)
        _pos("ruff", self.ruff)
        _pos("rufc", self.rufc)

        # Construct and populate the C++ object
        g       = _fred.FuelRodGeometry()
        g.nf    = self.nf
        g.nc    = self.nc
        g.nz    = self.nz
        g.rfi0  = self.rfi0
        g.rfo0  = self.rfo0
        g.rci0  = self.rci0
        g.rco0  = self.rco0
        g.dz0   = dz0
        g.vgp   = self.vgp
        g.ruff  = self.ruff
        g.rufc  = self.rufc
        g.build()
        self._cpp_obj = g
        return self

    @property
    def _cpp(self):
        """The underlying C++ fred.FuelRodGeometry (available after build())."""
        if self._cpp_obj is None:
            raise RuntimeError(
                "FuelRodGeometry.build() must be called before the geometry "
                "can be used in a solver"
            )
        return self._cpp_obj


# ---------------------------------------------------------------------------
# Material factories — validate constructor arguments, return C++ objects
# ---------------------------------------------------------------------------

def UO2(reference_density=10400.0):
    """
    UO2 fuel pellet.

    reference_density : as-fabricated density [kg/m³], default 10400 (≈95% TD)
    """
    _pos("reference_density", reference_density)
    return _fred.UO2(reference_density)


def MOX(pu_content, reference_density=-1.0):
    """
    MOX fuel pellet (Philipponneau 1992 thermal conductivity).

    pu_content        : Pu mole fraction, must be in (0.01, 0.50)
    reference_density : as-fabricated density [kg/m³]; ≤ 0 → 95 % of
                        theoretical density (auto)
    """
    _range_incl("pu_content", pu_content, 0.01, 0.50)
    if reference_density > 0:
        _pos("reference_density", reference_density)
    return _fred.MOX(pu_content, reference_density)


def AIM1(reference_density=7900.0):
    """AIM1 austenitic stainless steel cladding."""
    _pos("reference_density", reference_density)
    return _fred.AIM1(reference_density)


def T91(reference_density=7750.0):
    """T91 ferritic-martensitic steel cladding (no irradiation)."""
    _pos("reference_density", reference_density)
    return _fred.T91(reference_density)


# --- FRED-OX materials ---

def FredOxMOX(pu_content, rof0, sto0=1.97):
    """
    MOX fuel for FRED-OX with burnup-dependent thermal conductivity.

    pu_content : Pu mole fraction, must be in (0.01, 0.50)
    rof0       : as-fabricated density [kg/m³], must be > 0
    sto0       : initial O/M stoichiometry, must be in [1.80, 2.00]
    """
    _range_incl("pu_content", pu_content, 0.01, 0.50)
    _pos("rof0", rof0)
    _range_incl("sto0", sto0, 1.80, 2.00)
    return _fred.FredOxMOX(pu_content, rof0, sto0)


def FredOxAIM1(reference_density=7900.0):
    """AIM1 cladding for FRED-OX (Luzzi irradiation creep model)."""
    _pos("reference_density", reference_density)
    return _fred.FredOxAIM1(reference_density)


def FredOxT91(reference_density=7750.0):
    """T91 cladding for FRED-OX (zero irradiation creep rate)."""
    _pos("reference_density", reference_density)
    return _fred.FredOxT91(reference_density)


# --- Simple materials (no constructor parameters to validate) ---

def DummyFuelPellet():
    """Constant-property dummy fuel pellet."""
    return _fred.DummyFuelPellet()


def DummyCladding():
    """Constant-property dummy cladding."""
    return _fred.DummyCladding()


def DummyGapMaterial():
    """He gap material (Lanning-Hann + radiation)."""
    return _fred.DummyGapMaterial()


def He():
    """Pure helium gap material."""
    return _fred.He()


def HeKrXe():
    """He-Kr-Xe mixed-gas gap material (50 % He, 4 % Kr, 46 % Xe)."""
    return _fred.HeKrXe()


def FredOxGapMaterial():
    """Gap material for FRED-OX (He/Kr/Xe mixture + fission gas)."""
    return _fred.FredOxGapMaterial()


# ---------------------------------------------------------------------------
# Internal helper — unwrap geometry
# ---------------------------------------------------------------------------

def _unwrap_geom(geometry):
    """Return the C++ FuelRodGeometry from either a wrapper or a raw C++ object."""
    if isinstance(geometry, FuelRodGeometry):
        return geometry._cpp
    return geometry


# ---------------------------------------------------------------------------
# FredRodSolver — validated wrapper
# ---------------------------------------------------------------------------

class FredRodSolver:
    """
    Validated Python wrapper for fred.FredRodSolver.

    Accepts either a fred_input.FuelRodGeometry wrapper or a raw
    fred.FuelRodGeometry C++ object as the first argument.

    All boundary-condition setters are validated before being forwarded to
    the C++ solver.
    """

    def __init__(self, geometry, fuel, clad, gap):
        cpp_geom = _unwrap_geom(geometry)
        self._nz = (geometry.nz if isinstance(geometry, FuelRodGeometry)
                    else cpp_geom.nz)
        self._nf = cpp_geom.nf
        self._nc = cpp_geom.nc
        self._solver = _fred.FredRodSolver(cpp_geom, fuel, clad, gap)

    # -- boundary conditions -----------------------------------------------

    def set_power_density_history(self, times, qv_W_m3):
        """
        Time-history of volumetric power density [W/m³].

        times    : list of time points [s], strictly increasing, len >= 2
        qv_W_m3  : list of power densities [W/m³], non-negative, same length
        """
        _check_time_history(
            "set_power_density_history",
            times, qv_W_m3,
            _all_nonneg_vals,
        )
        self._solver.set_power_density_history(list(times), list(qv_W_m3))

    def set_coolant_temperature(self, times, T_K):
        """
        Time-history of outer cladding surface temperature [K].

        times : strictly increasing time points [s]
        T_K   : temperatures in Kelvin (must be positive)
        """
        _check_time_history(
            "set_coolant_temperature",
            times, T_K,
            _all_positive_K,
        )
        self._solver.set_coolant_temperature(list(times), list(T_K))

    def set_coolant_pressure(self, pcool_MPa):
        """Constant coolant pressure [MPa], must be positive."""
        _pos("pcool_MPa", pcool_MPa)
        self._solver.set_coolant_pressure(pcool_MPa)

    def set_coolant_pressure_history(self, times, pcool_MPa):
        """Time-varying coolant pressure history [MPa], all values must be positive."""
        _check_time_history(
            "set_coolant_pressure_history",
            times, pcool_MPa,
            _all_positive_MPa,
        )
        self._solver.set_coolant_pressure_history(list(times), list(pcool_MPa))

    def set_initial_temperature(self, T0_K):
        """Uniform initial temperature [K], must be positive."""
        _pos("T0_K", T0_K)
        self._solver.set_initial_temperature(T0_K)

    def set_initial_gas_pressure(self, gpres0_MPa):
        """Initial fill-gas pressure [MPa], must be positive."""
        _pos("gpres0_MPa", gpres0_MPa)
        self._solver.set_initial_gas_pressure(gpres0_MPa)

    def set_tolerances(self, rtol, atol):
        """IDA relative and absolute tolerances (both must be positive)."""
        _pos("rtol", rtol)
        _pos("atol", atol)
        self._solver.set_tolerances(rtol, atol)

    def set_enable_heat_conduction(self, enable):
        """Enable (True) or disable (False) the heat conduction block."""
        self._solver.set_enable_heat_conduction(bool(enable))

    def set_enable_stress_strain(self, enable):
        """Enable (True) or disable (False) the stress-strain block."""
        self._solver.set_enable_stress_strain(bool(enable))

    # -- run ---------------------------------------------------------------

    def run(self, tend, dtout, output_file=None, all_steps=False):
        """
        Run from t=0 to tend [s], writing output every dtout [s].

        tend        : positive end time [s]
        dtout       : positive output interval [s], must be <= tend
        output_file : if given, write results to this HDF5 file after the run
                      (e.g. 'results.h5').  Requires h5py.
        all_steps   : if True, record every internal IDA step in addition to
                      the dtout checkpoints.  Produces a much larger file and
                      is primarily for debugging.
        """
        _pos("tend", tend)
        _pos("dtout", dtout)
        if dtout > tend:
            raise ValueError(
                f"dtout ({dtout} s) must be <= tend ({tend} s)"
            )
        if output_file is not None:
            self._solver.set_output_file(str(output_file))
        self._solver.run(tend, dtout, bool(all_steps))

    # -- result accessors (delegated, unchanged) ---------------------------

    def time_points(self):            return self._solver.time_points()
    def temperatures(self):           return self._solver.temperatures()
    def peak_fuel_temperature(self):  return self._solver.peak_fuel_temperature()
    def clad_outer_hoop_stress(self): return self._solver.clad_outer_hoop_stress()
    def clad_outer_radius(self):      return self._solver.clad_outer_radius()
    def gap_width(self):              return self._solver.gap_width()
    def contact_pressure(self):       return self._solver.contact_pressure()
    def fuel_outer_radius(self):      return self._solver.fuel_outer_radius()
    def clad_inner_radius(self):      return self._solver.clad_inner_radius()
    def y_out(self):                  return self._solver.y_out()
    def yp_out(self):                 return self._solver.yp_out()


# ---------------------------------------------------------------------------
# FredOxSolver — validated wrapper
# ---------------------------------------------------------------------------

class FredOxSolver:
    """
    Validated Python wrapper for fred.FredOxSolver.

    Accepts either a fred_input.FuelRodGeometry wrapper or a raw
    fred.FuelRodGeometry C++ object as the first argument.
    The mox and gap arguments must be FredOxMOX / FredOxGapMaterial C++
    objects (as returned by the factory functions above).
    """

    def __init__(self, geometry, mox, clad, gap):
        cpp_geom = _unwrap_geom(geometry)
        self._nz = (geometry.nz if isinstance(geometry, FuelRodGeometry)
                    else cpp_geom.nz)
        self._nf = cpp_geom.nf
        self._nc = cpp_geom.nc
        self._solver = _fred.FredOxSolver(cpp_geom, mox, clad, gap)

    # -- power density -----------------------------------------------------

    def set_power_density_history(self, times, qqv_W_m3):
        """
        Uniform volumetric power density [W/m³] vs time (broadcast to all layers).

        times     : strictly increasing time points [s], len >= 2
        qqv_W_m3  : non-negative power densities [W/m³]
        """
        _check_time_history(
            "set_power_density_history",
            times, qqv_W_m3,
            _all_nonneg_vals,
        )
        self._solver.set_power_density_history(list(times), list(qqv_W_m3))

    def set_power_density_history_per_layer(self, times, qqv_per_layer):
        """
        Per-layer volumetric power density [W/m³] vs time.

        times          : strictly increasing time points [s], len >= 2
        qqv_per_layer  : list of nz lists; qqv_per_layer[j] is the power
                         history for axial layer j (len == len(times),
                         all values non-negative)
        """
        _check_per_layer_history(
            "set_power_density_history_per_layer",
            times, qqv_per_layer, self._nz,
        )
        for j, vals in enumerate(qqv_per_layer):
            _all_nonneg_vals(
                f"set_power_density_history_per_layer layer {j}", vals
            )
        self._solver.set_power_density_history_per_layer(
            list(times), [list(v) for v in qqv_per_layer]
        )

    # -- coolant temperature -----------------------------------------------

    def set_coolant_temperature(self, times, T_K):
        """Uniform outer cladding temperature [K] vs time (all layers)."""
        _check_time_history(
            "set_coolant_temperature",
            times, T_K,
            _all_positive_K,
        )
        self._solver.set_coolant_temperature(list(times), list(T_K))

    def set_coolant_temperature_per_layer(self, times, T_per_layer):
        """
        Per-layer outer cladding temperature [K] vs time.

        T_per_layer[j] : temperature history for axial layer j (all > 0 K)
        """
        _check_per_layer_history(
            "set_coolant_temperature_per_layer",
            times, T_per_layer, self._nz,
        )
        for j, vals in enumerate(T_per_layer):
            _all_positive_K(
                f"set_coolant_temperature_per_layer layer {j}", vals
            )
        self._solver.set_coolant_temperature_per_layer(
            list(times), [list(v) for v in T_per_layer]
        )

    # -- plenum temperature ------------------------------------------------

    def set_plenum_temperature_history(self, times, T_K):
        """Gas plenum temperature [K] vs time."""
        _check_time_history(
            "set_plenum_temperature_history",
            times, T_K,
            _all_positive_K,
        )
        self._solver.set_plenum_temperature_history(list(times), list(T_K))

    # -- pressure ----------------------------------------------------------

    def set_coolant_pressure(self, pcool_MPa):
        """Coolant pressure [MPa], must be positive."""
        _pos("pcool_MPa", pcool_MPa)
        self._solver.set_coolant_pressure(pcool_MPa)

    # -- initial conditions ------------------------------------------------

    def set_initial_temperature(self, T0_K):
        """Uniform initial temperature [K], must be positive."""
        _pos("T0_K", T0_K)
        self._solver.set_initial_temperature(T0_K)

    def set_initial_gas_pressure(self, gpres0_MPa):
        """Initial fill-gas pressure [MPa], must be positive."""
        _pos("gpres0_MPa", gpres0_MPa)
        self._solver.set_initial_gas_pressure(gpres0_MPa)

    # -- physics options ---------------------------------------------------

    def set_swelling_multiplier(self, fswelmlt):
        """
        Fuel swelling multiplier (1.0 = MATPRO full; typical: 0.8).
        Must be positive.
        """
        _pos("fswelmlt", fswelmlt)
        self._solver.set_swelling_multiplier(fswelmlt)

    def set_tolerances(self, rtol, atol):
        """IDA relative and absolute tolerances (both must be positive)."""
        _pos("rtol", rtol)
        _pos("atol", atol)
        self._solver.set_tolerances(rtol, atol)

    def set_enable_heat_conduction(self, enable):
        self._solver.set_enable_heat_conduction(bool(enable))

    def set_enable_stress_strain(self, enable):
        self._solver.set_enable_stress_strain(bool(enable))

    # -- run ---------------------------------------------------------------

    def run(self, tend, dtout, output_file=None, all_steps=False):
        """
        Run from t=0 to tend [s], writing output every dtout [s].

        tend        : positive end time [s]
        dtout       : positive output interval [s], must be <= tend
        output_file : if given, write results to this HDF5 file after the run
                      (e.g. 'results.h5').  Requires h5py.
        all_steps   : if True, record every internal IDA step in addition to
                      the dtout checkpoints.  Produces a much larger file and
                      is primarily for debugging.
        """
        _pos("tend", tend)
        _pos("dtout", dtout)
        if dtout > tend:
            raise ValueError(
                f"dtout ({dtout} s) must be <= tend ({tend} s)"
            )
        if output_file is not None:
            self._solver.set_output_file(str(output_file))
        self._solver.run(tend, dtout, bool(all_steps))

    # -- result accessors (delegated) --------------------------------------

    def time_points(self):           return self._solver.time_points()
    def temperatures(self):          return self._solver.temperatures()
    def peak_fuel_temperature(self): return self._solver.peak_fuel_temperature()
    def gas_pressure(self):          return self._solver.gas_pressure()
    def fg_generated(self):          return self._solver.fg_generated()
    def fg_released(self):           return self._solver.fg_released()
    def gap_width(self):             return self._solver.gap_width()
    def burnup(self):                return self._solver.burnup()
    def y_out(self):                 return self._solver.y_out()
    def yp_out(self):                return self._solver.yp_out()
