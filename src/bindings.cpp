#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_fred_rod  (py::module_& m);
void bind_fred_ox   (py::module_& m);
void bind_fred_m_na (py::module_& m);

PYBIND11_MODULE(fred, m) {
    m.doc() = "FRED 2.0 fuel performance platform — Python interface";
    // bind_fred_rod registers platform abstract types (FuelPelletMaterial,
    // CladdingMaterial, GapMaterial, FuelRodGeometry) first because OX and
    // M-Na concrete subclasses depend on them being already registered.
    bind_fred_rod  (m);
    bind_fred_ox   (m);
    bind_fred_m_na (m);
}
