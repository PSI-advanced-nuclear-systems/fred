#include <pybind11/pybind11.h>

namespace py = pybind11;
void bind_fred_rod  (py::module_& m);  // registers platform abstract types
void bind_fred_m_na (py::module_& m);

PYBIND11_MODULE(_fred_m_na, m) {
    m.doc() = "FRED-M-Na: U-Pu-Zr metallic fuel / HT-9 cladding / sodium-cooled fast reactor";
    bind_fred_rod(m);
    bind_fred_m_na(m);
}
