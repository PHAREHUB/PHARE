#ifndef PHARE_PY_MHD_HPP
#define PHARE_PY_MHD_HPP

#include <cstdint>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <utility>

#include "amr/solvers/time_integrator/euler_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk2_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk3_integrator.hpp"
#include "amr/solvers/time_integrator/ssprk4_5_integrator.hpp"

#include "core/numerics/reconstructions/constant.hpp"
#include "core/numerics/reconstructions/linear.hpp"
#include "core/numerics/reconstructions/weno3.hpp"
#include "core/numerics/reconstructions/wenoz.hpp"
#include "core/numerics/reconstructions/mp5.hpp"

#include "core/numerics/slope_limiters/min_mod.hpp"
#include "core/numerics/slope_limiters/van_leer.hpp"

#include "core/numerics/riemann_solvers/rusanov.hpp"
#include "core/numerics/riemann_solvers/hll.hpp"

#include "core/numerics/MHD_equations/MHD_equations.hpp"
#include "python3/mhd_defaults/default_mhd_registerer.hpp"

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

enum class TimeIntegratorType : uint8_t { Euler, TVDRK2, TVDRK3, SSPRK4_5, count };
enum class ReconstructionType : uint8_t { Constant, Linear, WENO3, WENOZ, MP5, count };
enum class SlopeLimiterType : uint8_t { VanLeer, MinMod, count };
enum class RiemannSolverType : uint8_t { Rusanov, HLL, count };

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
struct TimeIntegratorSelector<TimeIntegratorType::SSPRK4_5>
{
    template<template<typename> typename FVmethod, typename MHDModel>
    using type = SSPRK4_5Integrator<FVmethod, MHDModel>;
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

template<>
struct ReconstructionSelector<ReconstructionType::MP5>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = MP5Reconstruction<GridLayout, SlopeLimiter>;
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



template<typename Dimension, typename InterpOrder, typename NbRefinedPart, TimeIntegratorType TI,
         ReconstructionType RC, SlopeLimiterType SL, RiemannSolverType RS, bool Hall,
         bool Resistivity, bool HyperResistivity>
class RegistererSelector
{
    template<template<typename> typename FVmethod, typename MHDModel>
    using TimeIntegrator = typename TimeIntegratorSelector<TI>::template type<FVmethod, MHDModel>;

    template<typename GridLayout, typename SlopeLimiter>
    using Reconstruction =
        typename ReconstructionSelector<RC>::template type<GridLayout, SlopeLimiter>;

    template<typename GridLayout, bool HallFlag>
    using RiemannSolver = typename RiemannSolverSelector<RS>::template type<GridLayout, HallFlag>;

    using SlopeLimiter = typename SlopeLimiterSelector<RC, SL>::type;

    using Registerer_t = Registerer<Dimension, InterpOrder, NbRefinedPart, TimeIntegrator,
                                    Reconstruction, SlopeLimiter, RiemannSolver, MHDEquations, Hall,
                                    Resistivity, HyperResistivity>;

public:
    static constexpr void declare_etc(py::module& m, std::string const& type_string)
    {
        if constexpr (!unwanted_simulators_())
            Registerer_t::declare_etc(m, type_string);
    }

    static constexpr void declare_sim(py::module& m, std::string const& type_string)
    {
        if constexpr (!unwanted_simulators_())
            Registerer_t::declare_sim(m, type_string);
    }

private:
    static constexpr bool unwanted_simulators_()
    {
        bool constexpr is_hyper_nohall = HyperResistivity && !Hall;

        return is_hyper_nohall;
    }
};



template<typename Dimension, typename InterpOrder, typename NbRefinedParts>
constexpr void declare_all_mhd_params(py::module& m)
{
    DefaultMHDRegisterer<Dimension, InterpOrder, NbRefinedParts>::declare_defaults(m);

    std::string type_name = "_" + std::to_string(Dimension{}()) + "_"
                            + std::to_string(InterpOrder{}()) + "_"
                            + std::to_string(NbRefinedParts{}());

    std::string variant_name = "euler_constant_rusanov";
    std::string full_type    = type_name + "_" + variant_name;
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::Euler,
    //                    ReconstructionType::Constant, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_sim(m,
    //                    full_type);
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::Euler,
    //                    ReconstructionType::Constant, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_etc(m,
    //                    full_type);

    // variant_name = "euler_constant_rusanov_hall";
    // full_type    = type_name + "_" + variant_name;
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::Euler,
    //                    ReconstructionType::Constant, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, true, false, false>::declare_sim(m,
    //                    full_type);
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::Euler,
    //                    ReconstructionType::Constant, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, true, false, false>::declare_etc(m,
    //                    full_type);

    // variant_name = "ssprk4_5_wenoz_rusanov_hall";
    // full_type    = type_name + "_" + variant_name;
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::SSPRK4_5,
    //                    ReconstructionType::WENOZ, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, true, false, false>::declare_sim(m,
    //                    full_type);
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::SSPRK4_5,
    //                    ReconstructionType::WENOZ, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, true, false, false>::declare_etc(m,
    //                    full_type);


    // variant_name = "ssprk4_5_mp5_rusanov";
    // full_type    = type_name + "_" + variant_name;

    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::SSPRK4_5,
    //                    ReconstructionType::WENOZ, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_sim(m,
    //                    full_type);

    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::SSPRK4_5,
    //                    ReconstructionType::WENOZ, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_etc(m,
    //                    full_type);

    // variant_name = "tvdrk3_mp5_rusanov";
    // full_type    = type_name + "_" + variant_name;

    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,
    //                    ReconstructionType::WENOZ, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_sim(m,
    //                    full_type);

    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,
    //                    ReconstructionType::WENOZ, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_etc(m,
    //                    full_type);

    // variant_name = "tvdrk3_wenoz_rusanov_hall";
    // full_type    = type_name + "_" + variant_name;
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,
    //                    ReconstructionType::WENOZ, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, true, false, false>::declare_sim(m,
    //                    full_type);
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,
    //                    ReconstructionType::WENOZ, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, true, false, false>::declare_etc(m,
    //                    full_type);
    //
    // variant_name = "tvdrk2_linear_vanleer_rusanov";
    // full_type    = type_name + "_" + variant_name;
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK2,
    //                    ReconstructionType::Linear, SlopeLimiterType::VanLeer,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_sim(m,
    //                    full_type);
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK2,
    //                    ReconstructionType::Linear, SlopeLimiterType::VanLeer,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_etc(m,
    //                    full_type);

    variant_name = "tvdrk2_linear_vanleer_rusanov_hall";
    full_type    = type_name + "_" + variant_name;

    RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK2,
                       ReconstructionType::Linear, SlopeLimiterType::VanLeer,
                       RiemannSolverType::Rusanov, true, false, false>::declare_sim(m, full_type);

    RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK2,
                       ReconstructionType::Linear, SlopeLimiterType::VanLeer,
                       RiemannSolverType::Rusanov, true, false, false>::declare_etc(m, full_type);

    // variant_name = "tvdrk3_weno3_rusanov";
    // full_type    = type_name + "_" + variant_name;
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,
    //                    ReconstructionType::WENO3, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_sim(m,
    //                    full_type);
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,
    //                    ReconstructionType::WENO3, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, false, false, false>::declare_etc(m,
    //                    full_type);

    // variant_name = "tvdrk3_weno3_rusanov_hall";
    // full_type    = type_name + "_" + variant_name;
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,
    //                    ReconstructionType::WENO3, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, true, false, false>::declare_sim(m,
    //                    full_type);
    //
    // RegistererSelector<Dimension, InterpOrder, NbRefinedParts, TimeIntegratorType::TVDRK3,
    //                    ReconstructionType::WENO3, SlopeLimiterType::count,
    //                    RiemannSolverType::Rusanov, true, false, false>::declare_etc(m,
    //                    full_type);

    // auto constexpr ti_tuple   = make_enum_tuple<TimeIntegratorType>();
    // auto constexpr rc_tuple   = make_enum_tuple<ReconstructionType>();
    // auto constexpr sl_tuple   = make_enum_tuple<SlopeLimiterType>();
    // auto constexpr rs_tuple   = make_enum_tuple<RiemannSolverType>();
    // auto constexpr bool_tuple = std::make_tuple(false, true);
    //
    // auto constexpr ti_size   = std::tuple_size_v<std::decay_t<decltype(ti_tuple)>>;
    // auto constexpr rc_size   = std::tuple_size_v<std::decay_t<decltype(rc_tuple)>>;
    // auto constexpr sl_size   = std::tuple_size_v<std::decay_t<decltype(sl_tuple)>>;
    // auto constexpr rs_size   = std::tuple_size_v<std::decay_t<decltype(rs_tuple)>>;
    // auto constexpr bool_size = 2ull;
    //
    // for_N<ti_size>([&](auto i_ti) {
    //     auto constexpr ti = std::get<i_ti>(ti_tuple);
    //     for_N<rc_size>([&](auto i_rc) {
    //         auto constexpr rc = std::get<i_rc>(rc_tuple);
    //         for_N<rs_size>([&](auto i_rs) {
    //             auto constexpr rs = std::get<i_rs>(rs_tuple);
    //             for_N<bool_size>([&](auto i_hall) {
    //                 auto constexpr hall = std::get<i_hall>(bool_tuple);
    //                 for_N<bool_size>([&](auto i_res) {
    //                     auto constexpr res = std::get<i_res>(bool_tuple);
    //                     for_N<bool_size>([&](auto i_hyper) {
    //                         auto constexpr hyper_res = std::get<i_hyper>(bool_tuple);
    //
    //                         // Reconstructions using slope limiters
    //                         if constexpr (rc == ReconstructionType::Linear)
    //                         {
    //                             for_N<sl_size>([&](auto i_sl) {
    //                                 auto constexpr sl = get<i_sl>(sl_tuple);
    //                                 std::string variant_name
    //                                     = (ti == TimeIntegratorType::Euler    ? "euler"
    //                                        : ti == TimeIntegratorType::TVDRK2 ? "tvdrk2"
    //                                                                           : "tvdrk3")
    //                                       + std::string("_")
    //                                       + (rc == ReconstructionType::Constant ? "constant"
    //                                          : rc == ReconstructionType::Linear ? "linear"
    //                                          : rc == ReconstructionType::WENO3  ? "weno3"
    //                                                                             : "wenoz")
    //                                       + std::string("_")
    //                                       + (sl == SlopeLimiterType::VanLeer ? "vanleer" :
    //                                       "minmod")
    //
    //                                       + std::string("_")
    //                                       + (rs == RiemannSolverType::Rusanov ? "rusanov" :
    //                                       "hll")
    //                                       + (hall ? "_hall" : "") + (res ? "_res" : "")
    //                                       + (hyper_res ? "_hyperres" : "");
    //
    //                                 std::string full_type = type_name + "_" + variant_name;
    //
    //                                 RegistererSelector<Dimension, InterpOrder, NbRefinedParts,
    //                                 ti,
    //                                                    rc, sl, rs, hall, res,
    //                                                    hyper_res>::declare_sim(m, full_type);
    //                                 RegistererSelector<Dimension, InterpOrder, NbRefinedParts,
    //                                 ti,
    //                                                    rc, sl, rs, hall, res,
    //                                                    hyper_res>::declare_etc(m, full_type);
    //                             });
    //                         }
    //                         else
    //                         {
    //                             std::string variant_name
    //                                 = (ti == TimeIntegratorType::Euler    ? "euler"
    //                                    : ti == TimeIntegratorType::TVDRK2 ? "tvdrk2"
    //                                                                       : "tvdrk3")
    //                                   + std::string("_")
    //                                   + (rc == ReconstructionType::Constant ? "constant"
    //                                      : rc == ReconstructionType::WENO3  ? "weno3"
    //                                                                         : "wenoz")
    //                                   + std::string("_")
    //                                   + (rs == RiemannSolverType::Rusanov ? "rusanov" : "hll")
    //                                   + (hall ? "_hall" : "") + (res ? "_res" : "")
    //                                   + (hyper_res ? "_hyperres" : "");
    //
    //                             std::string full_type = type_name + "_" + variant_name;
    //
    //                             auto constexpr nosl = SlopeLimiterType::count; // returns void
    //
    //
    //                             RegistererSelector<Dimension, InterpOrder, NbRefinedParts, ti,
    //                             rc,
    //                                                nosl, rs, hall, res,
    //                                                hyper_res>::declare_sim(m, full_type);
    //
    //                             RegistererSelector<Dimension, InterpOrder, NbRefinedParts, ti,
    //                             rc,
    //                                                nosl, rs, hall, res,
    //                                                hyper_res>::declare_etc(m, full_type);
    //                         }
    //                     });
    //                 });
    //             });
    //         });
    //     });
    // });
}

} // namespace PHARE::pydata

#endif
