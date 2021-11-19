
#include "python3/pybind_def.h"

namespace py = pybind11;

namespace PHARE::pydata
{
PYBIND11_MODULE(cpp_etc, m)
{
    py::class_<core::Span<double>, std::shared_ptr<core::Span<double>>>(m, "Span");
    py::class_<PyArrayWrapper<double>, std::shared_ptr<PyArrayWrapper<double>>, core::Span<double>>(
        m, "PyWrapper");

    m.def("makePyArrayWrapper", makePyArrayWrapper<double>);
}
} // namespace PHARE::pydata
