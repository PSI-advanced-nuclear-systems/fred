"""
Shared geometry class for all FRED applications.
Not imported directly by users — re-exported from each app module.
"""


class FuelRodGeometry:
    """
    Fuel rod geometry specification.

    Set all attributes then call build() to validate.

    Required
    --------
    nf   : int >= 2          number of fuel radial nodes
    nc   : int >= 2          number of cladding radial nodes
    nz   : int >= 1          number of axial layers
    rfi0 : float >= 0        inner fuel radius [m] (0 for solid pellet)
    rfo0 : float > rfi0      outer fuel radius [m]
    rci0 : float > rfo0      inner cladding radius [m]
    rco0 : float > rci0      outer cladding radius [m]
    dz0  : list[float]       axial layer heights [m], length nz, all > 0

    Optional (have defaults)
    ------------------------
    vgp  : float >= 0        gas plenum volume [m³] (default 0.0)
    ruff : float > 0         fuel outer roughness [m] (default 1 µm)
    rufc : float > 0         cladding inner roughness [m] (default 1 µm)
    """

    def __init__(self):
        self.nf   = None
        self.nc   = None
        self.nz   = None
        self.rfi0 = None
        self.rfo0 = None
        self.rci0 = None
        self.rco0 = None
        self.dz0  = None
        self.vgp  = 0.0
        self.ruff = 1.0e-6
        self.rufc = 1.0e-6
        self._built = False
        self._cpp_cache = None

    def build(self):
        """Validate all geometry parameters. Must be called before passing to a solver."""
        for attr in ("nf", "nc", "nz", "rfi0", "rfo0", "rci0", "rco0", "dz0"):
            if getattr(self, attr) is None:
                raise ValueError(f"FuelRodGeometry.{attr} is not set")

        if not isinstance(self.nf, int) or self.nf < 2:
            raise ValueError(f"nf must be an integer >= 2, got {self.nf!r}")
        if not isinstance(self.nc, int) or self.nc < 2:
            raise ValueError(f"nc must be an integer >= 2, got {self.nc!r}")
        if not isinstance(self.nz, int) or self.nz < 1:
            raise ValueError(f"nz must be an integer >= 1, got {self.nz!r}")

        if self.rfi0 < 0:
            raise ValueError(f"rfi0 must be >= 0, got {self.rfi0!r}")
        if self.rfo0 <= self.rfi0:
            raise ValueError(
                f"rfo0 ({self.rfo0:.6g} m) must be > rfi0 ({self.rfi0:.6g} m)"
            )
        if self.rci0 <= self.rfo0:
            raise ValueError(
                f"rci0 ({self.rci0:.6g} m) must be > rfo0 ({self.rfo0:.6g} m) "
                "(as-fabricated gap must be positive)"
            )
        if self.rco0 <= self.rci0:
            raise ValueError(
                f"rco0 ({self.rco0:.6g} m) must be > rci0 ({self.rci0:.6g} m)"
            )

        dz0 = list(self.dz0)
        if len(dz0) != self.nz:
            raise ValueError(f"dz0 must have nz={self.nz} entries, got {len(dz0)}")
        for i, dz in enumerate(dz0):
            if dz <= 0:
                raise ValueError(f"dz0[{i}] must be positive [m], got {dz!r}")

        if self.vgp < 0:
            raise ValueError(f"vgp must be >= 0, got {self.vgp!r}")
        if self.ruff <= 0:
            raise ValueError(f"ruff must be > 0, got {self.ruff!r}")
        if self.rufc <= 0:
            raise ValueError(f"rufc must be > 0, got {self.rufc!r}")

        self._built = True
        self._cpp_cache = None  # invalidate any cached C++ object
        return self

    def _make_cpp(self, cpp_module):
        """
        Build and cache the C++ geometry object using the given pybind11 module.
        Called internally by each app solver — not part of the user API.
        """
        if not self._built:
            raise RuntimeError(
                "FuelRodGeometry.build() must be called before passing to a solver"
            )
        if self._cpp_cache is None:
            g = cpp_module.FuelRodGeometry()
            g.nf   = self.nf
            g.nc   = self.nc
            g.nz   = self.nz
            g.rfi0 = self.rfi0
            g.rfo0 = self.rfo0
            g.rci0 = self.rci0
            g.rco0 = self.rco0
            g.dz0  = list(self.dz0)
            g.vgp  = self.vgp
            g.ruff = self.ruff
            g.rufc = self.rufc
            g.build()
            self._cpp_cache = g
        return self._cpp_cache
