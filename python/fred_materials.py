"""
Abstract base classes for FRED material objects.
Concrete subclasses are defined in each app module (fred_rod.py, fred_ox.py,
fred_m_na.py) and re-exported at the app level.
"""


class FuelPellet:
    """Abstract base for fuel pellet materials."""

    def _make_cpp(self, cpp_module):
        raise NotImplementedError(f"{type(self).__name__}._make_cpp not implemented")


class Cladding:
    """Abstract base for cladding materials."""

    def _make_cpp(self, cpp_module):
        raise NotImplementedError(f"{type(self).__name__}._make_cpp not implemented")


class GapFill:
    """Abstract base for gap fill materials."""

    def _make_cpp(self, cpp_module):
        raise NotImplementedError(f"{type(self).__name__}._make_cpp not implemented")


def _validate_density(name, value):
    if value is not None and value <= 0:
        raise ValueError(f"{name} reference_density must be positive [kg/m³], got {value!r}")
