/*
 This file is for functions which are not coupled to the rest of the phare main simulator code

 ext = extraneous, could be changed to "etc"
*/

#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/chrono.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"

#include "core/utilities/types.hpp"
#include "python3/pybind_def.h"

namespace py = pybind11;

namespace PHARE::pydata
{
PYBIND11_MODULE(cpp_ext, m)
{
    py::class_<core::Span<double>, std::shared_ptr<core::Span<double>>>(m, "Span");
    py::class_<PyArrayWrapper<double>, std::shared_ptr<PyArrayWrapper<double>>, core::Span<double>>(
        m, "PyWrapper");

    m.def("makePyArrayWrapper", makePyArrayWrapper<double>);
}


} // namespace PHARE::pydata
