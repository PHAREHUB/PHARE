#ifndef PHARE_PYTHON_CPP_SIMULATOR_HPP
#define PHARE_PYTHON_CPP_SIMULATOR_HPP


#ifndef PHARE_SIM_STR
// Default template parameters for IDE/tooling parsing only
// This will be overridden at compile time with actual values (e.g., "1,1,2")
// via -DPHARE_SIM_STR=<dimension>,<interp_order>,<refined_particles>
#define PHARE_SIM_STR 1, 1, 2 // mostly for clangformat - errors in cpp file if define is missing
#endif

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

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


template<typename Simulator, typename PyClass>
void declareSimulator(PyClass&& sim)
{
    // Bind common simulator methods available to all template instantiations
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

template<typename Sim>
void inline declare_etc(py::module& m)
{
    // Extract template parameters from PHARE_SIM_STR macro (e.g., "1,1,2")
    constexpr auto opts = SimOpts{PHARE_SIM_STR};

    // Declare DataWrangler for accessing simulation data
    using DW         = DataWrangler<opts>;
    std::string name = "DataWrangler";

    py::class_<DW, py::smart_holder>(m, name.c_str())
        .def(py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def(py::init<std::shared_ptr<ISimulator> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def("sync_merge", &DW::sync_merge)
        .def("getPatchLevel", &DW::getPatchLevel)
        .def("getNumberOfLevels", &DW::getNumberOfLevels);

    // Declare PatchLevel for accessing patch-level data (fields, particles, etc.)
    using PL = PatchLevel<opts>;
    name     = "PatchLevel";
    py::class_<PL, py::smart_holder>(m, name.c_str())
        .def("getEM", &PL::getEM)
        .def("getE", &PL::getE)
        .def("getB", &PL::getB)
        .def("getBx", &PL::getBx)
        .def("getBy", &PL::getBy)
        .def("getBz", &PL::getBz)
        .def("getEx", &PL::getEx)
        .def("getEy", &PL::getEy)
        .def("getEz", &PL::getEz)
        .def("getVix", &PL::getVix)
        .def("getViy", &PL::getViy)
        .def("getViz", &PL::getViz)
        .def("getDensity", &PL::getDensity)
        .def("getBulkVelocity", &PL::getBulkVelocity)
        .def("getPopDensities", &PL::getPopDensities)
        .def("getPopFluxes", &PL::getPopFlux)
        .def("getFx", &PL::getFx)
        .def("getFy", &PL::getFy)
        .def("getFz", &PL::getFz)
        .def("getParticles", &PL::getParticles, py::arg("userPopName") = "all");

    // Declare Splitter for particle refinement operations
    using _Splitter
        = PHARE::amr::Splitter<core::DimConst<Sim::dimension>, core::InterpConst<Sim::interp_order>,
                               core::RefinedParticlesConst<Sim::nbRefinedPart>>;
    name = "Splitter";

    py::class_<_Splitter, py::smart_holder>(m, name.c_str())
        .def(py::init<>())
        .def_property_readonly_static("weight", [](py::object) { return _Splitter::weight; })
        .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });

    name = "split_pyarray_particles";
    m.def(name.c_str(), splitPyArrayParticles<_Splitter>);
}

void inline declare_macro_sim(py::module& m)
{
    // Instantiate Simulator with template parameters from PHARE_SIM_STR macro
    using Sim = Simulator<SimOpts{PHARE_SIM_STR}>;

    std::string name = "Simulator";
    declareSimulator<Sim>(
        py::class_<Sim, py::smart_holder>(m, name.c_str())
            // Expose template parameters as static properties for introspection
            .def_property_readonly_static("dims", [](py::object) { return Sim::dimension; })
            .def_property_readonly_static("interp_order",
                                          [](py::object) { return Sim::interp_order; })
            .def_property_readonly_static("refined_particle_nbr",
                                          [](py::object) { return Sim::nbRefinedPart; }));

    name = "make_simulator";
    m.def(name.c_str(), [](std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
        return makeSimulator<Sim>(hier);
    });


    declare_etc<Sim>(m);
}


} // namespace PHARE::pydata


// https://stackoverflow.com/a/51061314/795574
// ASAN detects leaks by default, even in system/third party libraries
inline char const* __asan_default_options()
{
    return "detect_leaks=0";
}



#endif /*PHARE_PYTHON_CPP_SIMULATOR_H*/
