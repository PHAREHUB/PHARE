#ifndef PHARE_PYTHON_CPP_SIMULATOR_HPP
#define PHARE_PYTHON_CPP_SIMULATOR_HPP


#include "phare/phare.hpp"

#include "core/utilities/mpi_utils.hpp"

#include "amr/wrappers/hierarchy.hpp"

#include "python3/particles.hpp"

#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/chrono.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"

#include <cstddef>

namespace py = pybind11;

namespace PHARE::pydata
{


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
}

template<typename Simulator, typename PyClass>
void declareSimulator(PyClass&& sim)
{
    sim.def("initialize", &Simulator::initialize)
        .def("advance", &Simulator::advance)
        .def("startTime", &Simulator::startTime)
        .def("currentTime", &Simulator::currentTime)
        .def("endTime", &Simulator::endTime)
        .def("timeStep", &Simulator::timeStep)
        .def("to_str", &Simulator::to_str)
        .def("domain_box", &Simulator::domainBox)
        .def("cell_width", &Simulator::cellWidth)
        .def("dump", &Simulator::dump, py::arg("timestamp"), py::arg("timestep"));
}

template<typename _dim, typename _interp, typename _nbRefinedPart>
void declare_etc(py::module& m)
{
    constexpr auto dim           = _dim{}();
    constexpr auto interp        = _interp{}();
    constexpr auto nbRefinedPart = _nbRefinedPart{}();

    std::string const type_string = "_" + std::to_string(dim) + "_" + std::to_string(interp) + "_"
                                    + std::to_string(nbRefinedPart);

    using _Splitter
        = PHARE::amr::Splitter<_dim, _interp, core::RefinedParticlesConst<nbRefinedPart>>;
    std::string name = "Splitter" + type_string;
    py::class_<_Splitter, std::shared_ptr<_Splitter>>(m, name.c_str())
        .def(py::init<>())
        .def_property_readonly_static("weight", [](py::object) { return _Splitter::weight; })
        .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });

    name = "split_pyarray_particles" + type_string;
    m.def(name.c_str(), splitPyArrayParticles<_Splitter>);
}

template<typename _dim, typename _interp, typename _nbRefinedPart>
void declare_sim(py::module& m)
{
    constexpr auto dim           = _dim{}();
    constexpr auto interp        = _interp{}();
    constexpr auto nbRefinedPart = _nbRefinedPart{}();

    std::string type_string = "_" + std::to_string(dim) + "_" + std::to_string(interp) + "_"
                              + std::to_string(nbRefinedPart);

    using Sim        = Simulator<dim, interp, nbRefinedPart>;
    std::string name = "Simulator" + type_string;
    declareSimulator<Sim>(
        py::class_<Sim, std::shared_ptr<Sim>>(m, name.c_str())
            .def_property_readonly_static("dims", [](py::object) { return Sim::dimension; })
            .def_property_readonly_static("interp_order",
                                          [](py::object) { return Sim::interp_order; })
            .def_property_readonly_static("refined_particle_nbr",
                                          [](py::object) { return Sim::nbRefinedPart; }));

    name = "make_simulator" + type_string;
    m.def(name.c_str(), [](std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
        return std::shared_ptr<Sim>{std::move(makeSimulator<dim, interp, nbRefinedPart>(hier))};
    });
}

template<typename Dimension, typename InterpOrder, typename... NbRefinedParts>
void declare_all(py::module& m, std::tuple<Dimension, InterpOrder, NbRefinedParts...> const&)
{
    core::apply(std::tuple<NbRefinedParts...>{}, [&](auto& nbRefinedPart) {
        declare_sim<Dimension, InterpOrder, std::decay_t<decltype(nbRefinedPart)>>(m);
        declare_etc<Dimension, InterpOrder, std::decay_t<decltype(nbRefinedPart)>>(m);
    });
}

void inline declare_essential(py::module& m)
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
