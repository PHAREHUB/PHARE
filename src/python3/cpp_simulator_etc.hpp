#ifndef PHARE_PYTHON_CPP_SIMULATOR_ETC_HPP
#define PHARE_PYTHON_CPP_SIMULATOR_ETC_HPP

#include "python3/cpp_simulator.hpp"

#include "python3/particles.hpp"     // IWYU pragma: keep
#include "python3/patch_level.hpp"   // IWYU pragma: keep
#include "python3/data_wrangler.hpp" // IWYU pragma: keep

namespace PHARE::pydata
{

template<auto opts>
void inline declare_etc_hybrid(py::module& m)
{
    using Sim = Simulator<opts>;

    using DW         = DataWrangler<opts>;
    std::string name = "DataWrangler";

    py::class_<DW, py::smart_holder>(m, name.c_str(), py::module_local())
        .def(py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy> const&>())
        .def("sync_merge", &DW::sync_merge)
        .def("getPatchLevel", &DW::getPatchLevel)
        .def("getNumberOfLevels", &DW::getNumberOfLevels);

    using PL = PatchLevel<opts>;
    name     = "PatchLevel";
    py::class_<PL, py::smart_holder>(m, name.c_str(), py::module_local())
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

// opts is a template parameter, so the discarded branch is never instantiated and an MHD-only
// build never processes the hybrid-only bindings.
template<auto opts>
void declare_etc(py::module& m)
{
    if constexpr (has_hybrid_v<opts>)
        declare_etc_hybrid<opts>(m);
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_CPP_SIMULATOR_ETC_HPP*/
