#ifndef PHARE_PYTHON_CPP_SIMULATOR_HPP
#define PHARE_PYTHON_CPP_SIMULATOR_HPP


#include "core/data/particles/particle_array_def.hpp"
#ifndef PHARE_SIM_STR
// mostly for clangformat - errors in cpp file if define is missing
#define PHARE_SIM_STR 1ull, 1ull, 2ull
#endif


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/def/phare_config.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"

#include "amr/samrai.hpp" // IWYU pragma: keep
#include "amr/wrappers/hierarchy.hpp"

#include "simulator/simulator.hpp" // IWYU pragma: keep

#include "python3/pybind_def.hpp" // IWYU pragma: keep
#include "simulator/simulator.hpp"

#include "pybind11/stl.h"        // IWYU pragma: keep
#include "pybind11/numpy.h"      // IWYU pragma: keep
#include "pybind11/chrono.h"     // IWYU pragma: keep
#include "pybind11/complex.h"    // IWYU pragma: keep
#include "pybind11/functional.h" // IWYU pragma: keep

#include "python3/particles.hpp"     // IWYU pragma: keep
#include "python3/patch_level.hpp"   // IWYU pragma: keep
#include "python3/data_wrangler.hpp" // IWYU pragma: keep
#include "python3/pybind_def.hpp"    // IWYU pragma: keep

#include "python3/patch_data.hpp"

#include "magic_enum/magic_enum_utility.hpp"

#include <cstddef>


namespace py = pybind11;

namespace PHARE::pydata
{
auto static constexpr resolve_simulator_options()
{
    using namespace PHARE;
    using namespace PHARE::core;
    using namespace PHARE::MHDOpts;
    return SimOpts{PHARE_SIM_STR};
}



template<typename Type, std::size_t dimension>
void declarePatchData(py::module& m, std::string key)
{
    using PatchDataType = PatchData<Type, dimension>;
    py::class_<PatchDataType>(m, key.c_str(), py::module_local())
        .def_readonly("patchID", &PatchDataType::patchID)
        .def_readonly("origin", &PatchDataType::origin)
        .def_readonly("lower", &PatchDataType::lower)
        .def_readonly("upper", &PatchDataType::upper)
        .def_readonly("nGhosts", &PatchDataType::nGhosts)
        .def_readonly("data", &PatchDataType::data);
}

template<SimOpts opts>
void declareParticles(py::module& m)
{
    using ParticleArray = SimOpts::Particles<opts>::value_type;
    std::string name    = "ParticleArray";
    py::class_<ParticleArray, std::shared_ptr<ParticleArray>>(m, name.c_str(), py::module_local());
    // .def_readonly("size", &ParticleArray::size);
    // .def("__getitem__",
    //      [](ParticleArray& self, std::size_t const idx) -> auto& { return self[idx]; })

    declarePatchData<ParticleArray, opts.dimension>(m, "PatchDataParticleArray");
    declarePatchData<ParticleArray*, opts.dimension>(m, "PatchDataParticleArrayPtr");

    declarePatchData<py_array_t<double>, opts.dimension>(m, "PatchPyArrayDouble");
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

template<typename Sim>
void inline declare_etc(py::module& m)
{
    constexpr auto opts = resolve_simulator_options();

    using DW         = DataWrangler<opts>;
    std::string name = "DataWrangler";

    py::class_<DW, py::smart_holder>(m, name.c_str())
        .def(py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def(py::init<std::shared_ptr<ISimulator> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def("sync", &DW::sync)
        .def("getMHDPatchLevel", &DW::getMHDPatchLevel)
        .def("getHybridPatchLevel", &DW::getHybridPatchLevel)
        .def("getNumberOfLevels", &DW::getNumberOfLevels);

    using HybPL = PatchLevel<typename Sim::HybridModel>;
    name        = "HybridPatchLevel";

    py::class_<HybPL, py::smart_holder>(m, name.c_str(), py::module_local())
        .def("getB", &HybPL::getB, py::arg("component"))
        .def("getE", &HybPL::getE, py::arg("component"))
        .def("getVi", &HybPL::getVi, py::arg("component"))
        .def("getN", &HybPL::getN, py::arg("pop_name"))
        .def("getNi", &HybPL::getNi)
        .def("getFlux", &HybPL::getFlux, py::arg("component"), py::arg("pop_name"))
        .def("getParticles", &HybPL::getParticles, py::arg("pop_name"));

    using MHDPL = PatchLevel<typename Sim::MHDModel>;
    name        = "MHDPatchLevel";
    py::class_<MHDPL, py::smart_holder>(m, name.c_str(), py::module_local());

    using _Splitter
        = PHARE::amr::Splitter<core::DimConst<Sim::dimension>, core::InterpConst<Sim::interp_order>,
                               core::RefinedParticlesConst<Sim::nbRefinedPart>>;
    name = "Splitter";

    py::class_<_Splitter, py::smart_holder>(m, name.c_str(), py::module_local())
        .def(py::init<>())
        .def_property_readonly_static("weight", [](py::object) { return _Splitter::weight; })
        .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });

    name = "split_pyarray_particles";
    m.def(name.c_str(), splitPyArrayParticles<_Splitter>);
}


void inline declare_macro_sim(py::module& m)
{
    constexpr auto opts = resolve_simulator_options();
    using Sim           = Simulator<opts>;

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



    declare_etc<Sim>(m);
    declareParticles<opts>(m);
}


} // namespace PHARE::pydata


// https://stackoverflow.com/a/51061314/795574
// ASAN detects leaks by default, even in system/third party libraries
inline char const* __asan_default_options()
{
    return "detect_leaks=0";
}



#endif /*PHARE_PYTHON_CPP_SIMULATOR_H*/
