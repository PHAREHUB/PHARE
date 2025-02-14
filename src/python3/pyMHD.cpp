#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "core/numerics/reconstructions/wenoz.hpp"
#include "core/utilities/types.hpp"
#include "tests/core/numerics/mock_mhd_simulator/test_mhd_solver.hpp"

#include "core/numerics/time_integrator/euler_integrator.hpp"
#include "core/numerics/time_integrator/tvdrk2_integrator.hpp"
#include "core/numerics/time_integrator/tvdrk3_integrator.hpp"

#include "core/numerics/reconstructions/constant.hpp"
#include "core/numerics/reconstructions/linear.hpp"
#include "core/numerics/reconstructions/weno3.hpp"

#include "core/numerics/slope_limiters/min_mod.hpp"
#include "core/numerics/slope_limiters/van_leer.hpp"

#include "core/numerics/riemann_solvers/rusanov.hpp"
#include "core/numerics/riemann_solvers/hll.hpp"

#include "core/numerics/MHD_equations/MHD_equations.hpp"


template<typename E>
using enum_underlying_type_t = std::underlying_type_t<E>;

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
using namespace PHARE::core;

template<std::size_t Constant>
using DimConst = PHARE::core::DimConst<Constant>;

template<std::size_t Constant>
using InterpConst = PHARE::core::InterpConst<Constant>;

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

template<std::size_t dim, std::size_t ord, TimeIntegratorType TI, ReconstructionType RC,
         SlopeLimiterType SL, RiemannSolverType RS, bool Hall, bool Resistivity,
         bool HyperResistivity>
class Registerer
{
    template<template<typename> typename FVmethod, typename MHDModel>
    using TimeIntegrator = typename TimeIntegratorSelector<TI>::template type<FVmethod, MHDModel>;

    template<typename GridLayout, typename SlopeLimiter>
    using Reconstruction =
        typename ReconstructionSelector<RC>::template type<GridLayout, SlopeLimiter>;

    template<typename GridLayout, bool HallFlag>
    using RiemannSolver = typename RiemannSolverSelector<RS>::template type<GridLayout, HallFlag>;

    using Limiter                              = typename SlopeLimiterSelector<RC, SL>::type;
    static constexpr std::size_t effective_ord = HyperResistivity ? ord + 1 : ord;

    using Simulator_t
        = MHDMockSimulator<dim, effective_ord, TimeIntegrator, Reconstruction, Limiter,
                           RiemannSolver, MHDEquations, Hall, Resistivity, HyperResistivity>;

public:
    static auto registerVariant(const std::string& type, py::module& m)
    {
        if (unwanted_simulators_())
            return;

        std::string name = "MHDMockSimulator_" + type;

        auto cls = py::class_<Simulator_t, std::shared_ptr<Simulator_t>>(m, name.c_str())
                       .def("advance", &Simulator_t::advance, py::arg("filename"),
                            py::arg("dumpfrequency"));

        name = "make_mhd_mock_simulator_" + type;
        m.def(name.c_str(), []() { return std::shared_ptr<Simulator_t>{std::move(create_())}; });
    };

private:
    static auto create_()
    {
        return std::make_unique<Simulator_t>(
            PHARE::initializer::PHAREDictHandler::INSTANCE().dict());
    }

    static bool unwanted_simulators_()
    {
        bool const is_hyper_nohall = HyperResistivity && !Hall;
        return is_hyper_nohall;
    }
};


template<std::size_t dim, std::size_t ord>
void registerSimulatorVariants(py::module& m)
{
    std::string type_name = std::to_string(dim) + "_" + std::to_string(ord);

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
                                    Registerer<dim, ord, ti, rc, sl, rs, hall, res,
                                               hyper_res>::registerVariant(full_type, m);
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
                                auto constexpr nosl   = SlopeLimiterType::count; // returns void
                                Registerer<dim, ord, ti, rc, nosl, rs, hall, res,
                                           hyper_res>::registerVariant(full_type, m);
                            }
                        });
                    });
                });
            });
        });
    });
}

PYBIND11_MODULE(pyMHD, m)
{
    /*registerSimulatorVariants<1, 1>(m);*/
    /*registerSimulatorVariants<2, 1>(m);*/
    /*registerSimulatorVariants<3, 1>(m);*/

    /*Registerer<2, 1, TimeIntegratorType::TVDRK2, ReconstructionType::Linear,*/
    /*           SlopeLimiterType::VanLeer, RiemannSolverType::Rusanov, false, false,*/
    /*           false>::registerVariant("2_1_tvdrk2_linear_vanleer_rusanov", m);*/

    /*Registerer<2, 1, TimeIntegratorType::TVDRK2, ReconstructionType::Constant,*/
    /*           SlopeLimiterType::count, RiemannSolverType::Rusanov, false, false,*/
    /*           false>::registerVariant("2_1_tvdrk2_constant_rusanov", m);*/

    Registerer<2, 2, TimeIntegratorType::TVDRK3, ReconstructionType::WENOZ, SlopeLimiterType::count,
               RiemannSolverType::Rusanov, false, false,
               false>::registerVariant("2_2_tvdrk3_wenoz_rusanov", m);
}
