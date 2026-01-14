#ifndef PHARE_CPP_MHD_PYTHON_REGISTERER_HPP
#define PHARE_CPP_MHD_PYTHON_REGISTERER_HPP

#include "phare_simulator_options.hpp"
#include "simulator/simulator.hpp"
#include "python3/mhd_resolver.hpp"
#include "python3/particles.hpp"
#include "python3/data_wrangler.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace PHARE::pydata
{

namespace py = pybind11;

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

template<typename Dimension, typename InterpOrder, typename NbRefinedPart,
         template<template<typename> typename, typename> typename TimeIntegrator,
         template<typename, typename> typename Reconstruction, typename SlopeLimiter,
         template<typename, bool> typename RiemannSolver,
         template<bool, bool, bool> typename Equations, bool Hall, bool Resistivity,
         bool HyperResistivity>
class Registerer
{
    static constexpr auto dim           = Dimension{}();
    static constexpr auto interp        = InterpOrder{}();
    static constexpr auto nbRefinedPart = NbRefinedPart{}();

    template<typename Model>
    using MHDTimeStepper_t =
        typename MHDResolver<TimeIntegrator, Reconstruction, SlopeLimiter, RiemannSolver, Equations,
                             Hall, Resistivity, HyperResistivity>::template TimeIntegrator_t<Model>;

    static constexpr SimOpts<MHDTimeStepper_t> opts{dim, interp, nbRefinedPart};

    using Sim = Simulator<opts>;
    using DW  = DataWrangler<opts>;

public:
    static constexpr void declare_etc(py::module& m, std::string const& full_type)
    {
        std::string name = "DataWrangler" + full_type;

        py::class_<DW, py::smart_holder>(m, name.c_str())
            .def(py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy> const&>())
            .def(py::init<std::shared_ptr<ISimulator> const&,
                          std::shared_ptr<amr::Hierarchy> const&>())
            .def("sync_merge", &DW::sync_merge)
            .def("getPatchLevel", &DW::getPatchLevel)
            .def("getNumberOfLevels", &DW::getNumberOfLevels);

        using PL = PatchLevel<opts>;
        name     = "PatchLevel_" + full_type;

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
    }

    static constexpr void declare_sim(py::module& m, std::string const& full_type)
    {
        std::string name = "Simulator" + full_type;

        declareSimulator<Sim>(
            py::class_<Sim, py::smart_holder>(m, name.c_str())
                .def_property_readonly_static("dims", [](py::object) { return Sim::dimension; })
                .def_property_readonly_static("interp_order",
                                              [](py::object) { return Sim::interp_order; })
                .def_property_readonly_static("refined_particle_nbr",
                                              [](py::object) { return Sim::nbRefinedPart; }));

        name = "make_simulator" + full_type;
        m.def(name.c_str(), [](std::shared_ptr<PHARE::amr::Hierarchy> const& hier) {
            return std::shared_ptr<Sim>{std::move(makeSimulator<Sim>(hier))};
        });
    }
};

} // namespace PHARE::pydata

#endif
