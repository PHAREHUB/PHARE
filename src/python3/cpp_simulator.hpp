#ifndef PHARE_PYTHON_CPP_SIMULATOR_HPP
#define PHARE_PYTHON_CPP_SIMULATOR_HPP

#include <vector>
#include <cstddef>

#include "core/def/phare_mpi.hpp"

#include "core/utilities/mpi_utils.hpp"
#include "core/data/particles/particle.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "phare/phare.hpp"
#include "simulator/simulator.hpp"

#include "python3/pybind_def.hpp"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/chrono.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"

#include "python3/particles.hpp"
#include "python3/patch_data.hpp"
#include "python3/patch_level.hpp"
#include "python3/data_wrangler.hpp"
#include "python3/cpp_mhd_parameters.hpp"
#include "python3/cpp_mhd_python_registerer.hpp"


namespace PHARE::pydata
{
template<typename Type, std::size_t dimension>
void declarePatchData(py::module& m, std::string key)
{
    using PatchDataType = PatchData<Type, dimension>;
    py::class_<PatchDataType>(m, key.c_str())
        .def_readonly("patchID", &PatchDataType::patchID)
        .def_readonly("origin", &PatchDataType::origin)
        .def_readonly("lower", &PatchDataType::lower)
        .def_readonly("upper", &PatchDataType::upper)
        .def_readonly("nGhosts", &PatchDataType::nGhosts)
        .def_readonly("data", &PatchDataType::data);
}

template<std::size_t dim>
void declareDim(py::module& m)
{
    using CP         = core::ContiguousParticles<dim>;
    std::string name = "ContiguousParticles_" + std::to_string(dim);
    py::class_<CP, std::shared_ptr<CP>>(m, name.c_str())
        .def(py::init<std::size_t>())
        .def_readwrite("iCell", &CP::iCell)
        .def_readwrite("delta", &CP::delta)
        .def_readwrite("weight", &CP::weight)
        .def_readwrite("charge", &CP::charge)
        .def_readwrite("v", &CP::v)
        .def("size", &CP::size);

    name = "PatchData" + name;
    declarePatchData<CP, dim>(m, name.c_str());
}

template<typename Dimension, typename InterpOrder, typename... NbRefinedParts>
void declare_all(py::module& m, std::tuple<Dimension, InterpOrder, NbRefinedParts...> const&)
{
    core::apply(std::tuple<NbRefinedParts...>{}, [&](auto& nbRefinedPart) {
        declare_all_mhd_params<Dimension, InterpOrder, std::decay_t<decltype(nbRefinedPart)>>(m);
    });
}

void declare_essential(py::module& m)
{
    py::class_<SamraiLifeCycle, std::shared_ptr<SamraiLifeCycle>>(m, "SamraiLifeCycle")
        .def(py::init<>())
        .def("reset", &SamraiLifeCycle::reset);

    py::class_<PHARE::amr::Hierarchy, std::shared_ptr<PHARE::amr::Hierarchy>>(m, "AMRHierarchy");
    m.def("make_hierarchy", []() { return PHARE::amr::Hierarchy::make(); });

    m.def("mpi_size", []() { return core::mpi::size(); });
    m.def("mpi_rank", []() { return core::mpi::rank(); });
    m.def("mpi_barrier", []() { core::mpi::barrier(); });
}


} // namespace PHARE::pydata


// https://stackoverflow.com/a/51061314/795574
// ASAN detects leaks by default, even in system/third party libraries
inline char const* __asan_default_options()
{
    return "detect_leaks=0";
}



#endif /*PHARE_PYTHON_CPP_SIMULATOR_H*/
