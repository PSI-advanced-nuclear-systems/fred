#include <pybind11/pybind11.h>

namespace py = pybind11;
void bind_fred_rod(py::module_& m);  // registers platform abstract types
void bind_fred_ox (py::module_& m);

PYBIND11_MODULE(_fred_ox, m) {
    m.doc() = "FRED-OX: MOX fuel pin with fission gas release and irradiation effects";
    bind_fred_rod(m);
    bind_fred_ox(m);
}
