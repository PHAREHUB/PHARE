#ifndef PHARE_PYTHON_CPP_SIMULATOR_HPP
#define PHARE_PYTHON_CPP_SIMULATOR_HPP


#ifndef PHARE_SIM_STR
#define PHARE_SIM_STR 1, 1, 2 // mostly for clangformat - errors in cpp file if define is missing
#endif

#include "phare_mpi.hpp" // IWYU pragma: keep

#include "amr/samrai.hpp" // IWYU pragma: keep
#include "amr/wrappers/hierarchy.hpp"

#include "simulator/simulator.hpp" // IWYU pragma: keep

#include "python3/pybind_def.hpp" // IWYU pragma: keep
#include "pybind11/stl.h"         // IWYU pragma: keep
#include "pybind11/numpy.h"       // IWYU pragma: keep
#include "pybind11/chrono.h"      // IWYU pragma: keep
#include "pybind11/complex.h"     // IWYU pragma: keep
#include "pybind11/functional.h"  // IWYU pragma: keep

#include "python3/particles.hpp"     // IWYU pragma: keep
#include "python3/patch_level.hpp"   // IWYU pragma: keep
#include "python3/data_wrangler.hpp" // IWYU pragma: keep


namespace py = pybind11;

namespace PHARE::pydata
{


auto static constexpr resolve_simulator_options()
{
    using namespace PHARE::MHDOpts;
    return SimOpts{PHARE_SIM_STR};
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
        .def("dump_diagnostics", &Simulator::dump_diagnostics, py::arg("timestamp"),
             py::arg("timestep"))
        .def("dump_restarts", &Simulator::dump_restarts, py::arg("timestamp"), py::arg("timestep"));
}

// every build gets the DataWrangler, with the patch level getters of the models it has
template<auto opts>
void inline declare_data_wrangler(py::module& m)
{
    using Sim = Simulator<opts>;
    using DW  = DataWrangler<opts>;

    auto dw
        = py::class_<DW, py::smart_holder>(m, "DataWrangler", py::module_local())
              .def(py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy> const&>())
              .def(py::init<std::shared_ptr<ISimulator> const&,
                            std::shared_ptr<amr::Hierarchy> const&>())
              .def("sync", &DW::sync)
              .def("getNumberOfLevels", &DW::getNumberOfLevels);

    if constexpr (has_hybrid_v<opts>)
        dw.def("getHybridPatchLevel", &DW::getHybridPatchLevel);
    if constexpr (has_mhd_v<opts>)
        dw.def("getMHDPatchLevel", &DW::getMHDPatchLevel);
}

template<auto opts>
void inline declare_etc_hybrid(py::module& m)
{
    using HybPL      = PatchLevel<typename solver::PHARE_Types<opts>::Hybrid::Model_t>;
    std::string name = "HybridPatchLevel";
    py::class_<HybPL, py::smart_holder>(m, name.c_str(), py::module_local())
        .def("getB", &HybPL::getB, py::arg("component"))
        .def("getE", &HybPL::getE, py::arg("component"))
        .def("getVi", &HybPL::getVi, py::arg("component"))
        .def("getN", &HybPL::getN, py::arg("pop_name"))
        .def("getNi", &HybPL::getNi)
        .def("getFlux", &HybPL::getFlux, py::arg("component"), py::arg("pop_name"))
        .def("getParticles", &HybPL::getParticles, py::arg("pop_name"));

    using _Splitter
        = PHARE::amr::Splitter<core::DimConst<opts.dimension>, core::InterpConst<opts.interp_order>,
                               core::RefinedParticlesConst<opts.nbRefinedPart>>;
    name = "Splitter";

    py::class_<_Splitter, py::smart_holder>(m, name.c_str(), py::module_local())
        .def(py::init<>())
        .def_property_readonly_static("weight", [](py::object) { return _Splitter::weight; })
        .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });

    name = "split_pyarray_particles";
    m.def(name.c_str(), splitPyArrayParticles<_Splitter>);
}

template<auto opts>
void inline declare_etc_mhd(py::module& m)
{
    using MHDPL      = PatchLevel<typename solver::PHARE_Types<opts>::MHD::Model_t>;
    std::string name = "MHDPatchLevel";
    py::class_<MHDPL, py::smart_holder>(m, name.c_str(), py::module_local())
        .def("getRho", &MHDPL::getRho)
        .def("getP", &MHDPL::getP)
        .def("getEtot", &MHDPL::getEtot)
        .def("getV", &MHDPL::getV, py::arg("component"))
        .def("getB", &MHDPL::getB, py::arg("component"))
        .def("getRhoV", &MHDPL::getRhoV, py::arg("component"))
        .def("getE", &MHDPL::getE, py::arg("component"))
        .def("getJ", &MHDPL::getJ, py::arg("component"));
}

// opts is a template parameter, so the discarded branches are never instantiated: a build only
// processes the bindings of the models it has.
template<auto opts>
void inline declare_etc(py::module& m)
{
    declare_data_wrangler<opts>(m);
    if constexpr (has_hybrid_v<opts>)
        declare_etc_hybrid<opts>(m);
    if constexpr (has_mhd_v<opts>)
        declare_etc_mhd<opts>(m);
}


void inline declare_macro_sim(py::module& m)
{
    using Sim = Simulator<resolve_simulator_options()>;

    std::string name = "Simulator";
    declareSimulator<Sim>(
        py::class_<Sim, py::smart_holder>(m, name.c_str())
            .def_property_readonly_static("dims", [](py::object) { return Sim::dimension; })
            .def_property_readonly_static("interp_order",
                                          [](py::object) { return Sim::interp_order; })
            .def_property_readonly_static("refined_particle_nbr",
                                          [](py::object) { return Sim::nbRefinedPart; }));

    name = "make_simulator";
    m.def(name.c_str(), [](std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
        return makeSimulator<Sim>(hier);
    });


    declare_etc<Sim::options>(m);
}


} // namespace PHARE::pydata


// https://stackoverflow.com/a/51061314/795574
// ASAN detects leaks by default, even in system/third party libraries
inline char const* __asan_default_options()
{
    return "detect_leaks=0";
}



#endif /*PHARE_PYTHON_CPP_SIMULATOR_H*/
