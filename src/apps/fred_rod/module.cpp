#include <pybind11/pybind11.h>

namespace py = pybind11;
void bind_fred_rod(py::module_& m);

PYBIND11_MODULE(_fred_rod, m) {
    m.doc() = "FRED-ROD: generic fuel pin — thermal conduction + stress-strain solver";
    bind_fred_rod(m);
}
