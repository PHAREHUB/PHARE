#ifndef PHARE_PY_MHD_HPP
#define PHARE_PY_MHD_HPP

#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <utility>

#include "python3/data_wrangler.hpp"
#include "python3/MHDResolver.hpp"
#include "python3/particles.hpp"

#include "core/utilities/types.hpp"
#include "simulator/simulator.hpp"

#include "amr/solvers/time_integrator/euler_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk2_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk3_integrator.hpp"

#include "core/numerics/reconstructions/constant.hpp"
#include "core/numerics/reconstructions/linear.hpp"
#include "core/numerics/reconstructions/weno3.hpp"
#include "core/numerics/reconstructions/wenoz.hpp"

#include "core/numerics/slope_limiters/min_mod.hpp"
#include "core/numerics/slope_limiters/van_leer.hpp"

#include "core/numerics/riemann_solvers/rusanov.hpp"
#include "core/numerics/riemann_solvers/hll.hpp"

#include "core/numerics/godunov_fluxes/godunov_fluxes.hpp"

#include "core/numerics/MHD_equations/MHD_equations.hpp"

namespace PHARE::pydata
{
template<typename E>
constexpr std::size_t enum_size()
{
    return static_cast<std::size_t>(E::count);
}

template<typename E, std::size_t... I>
constexpr auto make_enum_tuple_impl(std::index_sequence<I...>)
{
    return std::make_tuple(static_cast<E>(I)...);
}

template<typename E>
constexpr auto make_enum_tuple()
{
    return make_enum_tuple_impl<E>(std::make_index_sequence<enum_size<E>()>{});
}

namespace py = pybind11;
using namespace core;
using namespace solver;

enum class TimeIntegratorType { Euler, TVDRK2, TVDRK3, count };
enum class ReconstructionType { Constant, Linear, WENO3, WENOZ, count };
enum class SlopeLimiterType { VanLeer, MinMod, count };
enum class RiemannSolverType { Rusanov, HLL, count };

template<TimeIntegratorType T>
struct TimeIntegratorSelector;

template<ReconstructionType T>
struct ReconstructionSelector;

template<ReconstructionType R, SlopeLimiterType S>
struct SlopeLimiterSelector;

template<RiemannSolverType T>
struct RiemannSolverSelector;

template<>
struct TimeIntegratorSelector<TimeIntegratorType::Euler>
{
    template<template<typename> typename FVmethod, typename MHDModel>
    using type = EulerIntegrator<FVmethod, MHDModel>;
};

template<>
struct TimeIntegratorSelector<TimeIntegratorType::TVDRK2>
{
    template<template<typename> typename FVmethod, typename MHDModel>
    using type = TVDRK2Integrator<FVmethod, MHDModel>;
};

template<>
struct TimeIntegratorSelector<TimeIntegratorType::TVDRK3>
{
    template<template<typename> typename FVmethod, typename MHDModel>
    using type = TVDRK3Integrator<FVmethod, MHDModel>;
};

template<>
struct ReconstructionSelector<ReconstructionType::Constant>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = ConstantReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<ReconstructionType::Linear>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = LinearReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<ReconstructionType::WENO3>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = WENO3Reconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<ReconstructionType::WENOZ>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = WENOZReconstruction<GridLayout, SlopeLimiter>;
};

template<ReconstructionType R, SlopeLimiterType S>
struct SlopeLimiterSelector
{
    using type = void;
};

template<>
struct SlopeLimiterSelector<ReconstructionType::Linear, SlopeLimiterType::VanLeer>
{
    using type = VanLeerLimiter;
};

template<>
struct SlopeLimiterSelector<ReconstructionType::Linear, SlopeLimiterType::MinMod>
{
    using type = MinModLimiter;
};

template<>
struct RiemannSolverSelector<RiemannSolverType::Rusanov>
{
    template<typename GridLayout, bool Hall>
    using type = Rusanov<GridLayout, Hall>;
};

template<>
struct RiemannSolverSelector<RiemannSolverType::HLL>
{
    template<typename GridLayout, bool Hall>
    using type = HLL<GridLayout, Hall>;
};


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

template<typename Dimension, typename InterpOrder, typename NbRefinedPart, TimeIntegratorType TI,
         ReconstructionType RC, SlopeLimiterType SL, RiemannSolverType RS, bool Hall,
         bool Resistivity, bool HyperResistivity>
class Registerer
{
    template<template<typename> typename FVmethod, typename MHDModel>
    using TimeIntegrator = typename TimeIntegratorSelector<TI>::template type<FVmethod, MHDModel>;

    template<typename GridLayout, typename SlopeLimiter>
    using Reconstruction =
        typename ReconstructionSelector<RC>::template type<GridLayout, SlopeLimiter>;

    template<typename GridLayout, bool HallFlag>
    using RiemannSolver = typename RiemannSolverSelector<RS>::template type<GridLayout, HallFlag>;

    using SlopeLimiter = typename SlopeLimiterSelector<RC, SL>::type;

    template<typename Model>
    using MHDTimeStepper_t =
        typename MHDResolver<TimeIntegrator, Reconstruction, SlopeLimiter, RiemannSolver,
                             MHDEquations, Hall, Resistivity,
                             HyperResistivity>::template TimeIntegrator_t<Model>;

    static constexpr auto dim           = Dimension{}();
    static constexpr auto interp        = InterpOrder{}();
    static constexpr auto nbRefinedPart = NbRefinedPart{}();

    using Sim = Simulator<dim, interp, nbRefinedPart, MHDTimeStepper_t>;
    using DW  = DataWrangler<dim, interp, nbRefinedPart, MHDTimeStepper_t>;

public:
    static constexpr void declare_etc(py::module& m, std::string const& type_string)
    {
        if constexpr (unwanted_simulators_())
            return;
        else
        {
            std::string name = "DataWrangler" + type_string;

            py::class_<DW, std::shared_ptr<DW>>(m, name.c_str())
                .def(
                    py::init<std::shared_ptr<Sim> const&, std::shared_ptr<amr::Hierarchy> const&>())
                .def(py::init<std::shared_ptr<ISimulator> const&,
                              std::shared_ptr<amr::Hierarchy> const&>())
                .def("sync_merge", &DW::sync_merge)
                .def("getPatchLevel", &DW::getPatchLevel)
                .def("getNumberOfLevels", &DW::getNumberOfLevels);

            using PL = PatchLevel<dim, interp, nbRefinedPart, MHDTimeStepper_t>;
            name     = "PatchLevel_" + type_string;

            py::class_<PL, std::shared_ptr<PL>>(m, name.c_str())
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

            using _Splitter = PHARE::amr::Splitter<Dimension, InterpOrder,
                                                   core::RefinedParticlesConst<nbRefinedPart>>;
            name            = "Splitter" + type_string;
            py::class_<_Splitter, std::shared_ptr<_Splitter>>(m, name.c_str())
                .def(py::init<>())
                .def_property_readonly_static("weight",
                                              [](py::object) { return _Splitter::weight; })
                .def_property_readonly_static("delta", [](py::object) { return _Splitter::delta; });

            name = "split_pyarray_particles" + type_string;
            m.def(name.c_str(), splitPyArrayParticles<_Splitter>);
        }
    }

    static constexpr void declare_sim(py::module& m, std::string const& type_string)
    {
        if constexpr (unwanted_simulators_())
            return;
        else
        {
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
                return std::shared_ptr<Sim>{
                    std::move(makeSimulator<dim, interp, nbRefinedPart, MHDTimeStepper_t>(hier))};
            });
        }
    }

private:
    static constexpr bool unwanted_simulators_()
    {
        /*bool constexpr is_hyper_nohall = HyperResistivity && !Hall;*/

        /*return is_hyper_nohall;*/

        bool constexpr compile_one
            = (TI == TimeIntegratorType::TVDRK3 && RC == ReconstructionType::WENOZ
               && SL == SlopeLimiterType::count && RS == RiemannSolverType::Rusanov
               && (Hall || !Hall) && !Resistivity && !HyperResistivity);

        return !compile_one;
    }
};


template<typename Dimension, typename InterpOrder, typename NbRefinedParts>
constexpr void declare_all_mhd_params(py::module& m)
{
    std::string type_name = "_" + std::to_string(Dimension{}()) + "_"
                            + std::to_string(InterpOrder{}()) + "_"
                            + std::to_string(NbRefinedParts{}());

    /*std::string variant_name = "tvdrk3_wenoz_rusanov";*/
    /*std::string full_type    = type_name + "_" + variant_name;*/
    /**/
    /*Registerer<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,*/
    /*           ReconstructionType::WENOZ, SlopeLimiterType::count, RiemannSolverType::Rusanov,*/
    /*           false, false, false>::declare_sim(m, full_type);*/
    /**/
    /*Registerer<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,*/
    /*           ReconstructionType::WENOZ, SlopeLimiterType::count, RiemannSolverType::Rusanov,*/
    /*           false, false, false>::declare_etc(m, full_type);*/

    auto constexpr ti_tuple   = make_enum_tuple<TimeIntegratorType>();
    auto constexpr rc_tuple   = make_enum_tuple<ReconstructionType>();
    auto constexpr sl_tuple   = make_enum_tuple<SlopeLimiterType>();
    auto constexpr rs_tuple   = make_enum_tuple<RiemannSolverType>();
    auto constexpr bool_tuple = std::make_tuple(false, true);

    auto constexpr ti_size   = std::tuple_size_v<std::decay_t<decltype(ti_tuple)>>;
    auto constexpr rc_size   = std::tuple_size_v<std::decay_t<decltype(rc_tuple)>>;
    auto constexpr sl_size   = std::tuple_size_v<std::decay_t<decltype(sl_tuple)>>;
    auto constexpr rs_size   = std::tuple_size_v<std::decay_t<decltype(rs_tuple)>>;
    auto constexpr bool_size = 2ull;

    for_N<ti_size>([&](auto i_ti) {
        auto constexpr ti = std::get<i_ti>(ti_tuple);
        for_N<rc_size>([&](auto i_rc) {
            auto constexpr rc = std::get<i_rc>(rc_tuple);
            for_N<rs_size>([&](auto i_rs) {
                auto constexpr rs = std::get<i_rs>(rs_tuple);
                for_N<bool_size>([&](auto i_hall) {
                    auto constexpr hall = std::get<i_hall>(bool_tuple);
                    for_N<bool_size>([&](auto i_res) {
                        auto constexpr res = std::get<i_res>(bool_tuple);
                        for_N<bool_size>([&](auto i_hyper) {
                            auto constexpr hyper_res = std::get<i_hyper>(bool_tuple);

                            // Reconstructions using slope limiters
                            if constexpr (rc == ReconstructionType::Linear)
                            {
                                for_N<sl_size>([&](auto i_sl) {
                                    auto constexpr sl = get<i_sl>(sl_tuple);
                                    std::string variant_name
                                        = (ti == TimeIntegratorType::Euler    ? "euler"
                                           : ti == TimeIntegratorType::TVDRK2 ? "tvdrk2"
                                                                              : "tvdrk3")
                                          + std::string("_")
                                          + (rc == ReconstructionType::Constant ? "constant"
                                             : rc == ReconstructionType::Linear ? "linear"
                                             : rc == ReconstructionType::WENO3  ? "weno3"
                                                                                : "wenoz")
                                          + std::string("_")
                                          + (sl == SlopeLimiterType::VanLeer ? "vanleer" : "minmod")

                                          + std::string("_")
                                          + (rs == RiemannSolverType::Rusanov ? "rusanov" : "hll")
                                          + (hall ? "_hall" : "") + (res ? "_res" : "")
                                          + (hyper_res ? "_hyperres" : "");

                                    std::string full_type = type_name + "_" + variant_name;

                                    Registerer<Dimension, InterpOrder, NbRefinedParts, ti, rc, sl,
                                               rs, hall, res, hyper_res>::declare_sim(m, full_type);

                                    Registerer<Dimension, InterpOrder, NbRefinedParts, ti, rc, sl,
                                               rs, hall, res, hyper_res>::declare_etc(m, full_type);
                                });
                            }
                            else
                            {
                                std::string variant_name
                                    = (ti == TimeIntegratorType::Euler    ? "euler"
                                       : ti == TimeIntegratorType::TVDRK2 ? "tvdrk2"
                                                                          : "tvdrk3")
                                      + std::string("_")
                                      + (rc == ReconstructionType::Constant ? "constant"
                                         : rc == ReconstructionType::WENO3  ? "weno3"
                                                                            : "wenoz")
                                      + std::string("_")
                                      + (rs == RiemannSolverType::Rusanov ? "rusanov" : "hll")
                                      + (hall ? "_hall" : "") + (res ? "_res" : "")
                                      + (hyper_res ? "_hyperres" : "");

                                std::string full_type = type_name + "_" + variant_name;

                                auto constexpr nosl = SlopeLimiterType::count; // returns void


                                Registerer<Dimension, InterpOrder, NbRefinedParts, ti, rc, nosl, rs,
                                           hall, res, hyper_res>::declare_sim(m, full_type);

                                Registerer<Dimension, InterpOrder, NbRefinedParts, ti, rc, nosl, rs,
                                           hall, res, hyper_res>::declare_etc(m, full_type);
                            }
                        });
                    });
                });
            });
        });
    });
}

} // namespace PHARE::pydata

#endif
