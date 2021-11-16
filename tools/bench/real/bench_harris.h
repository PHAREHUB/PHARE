#ifndef PHARE_BENCH_HARRIS_BENCH_H
#define PHARE_BENCH_HARRIS_BENCH_H

#include "mkn/kul/dbg.hpp"
#include "mkn/kul/log.hpp"

#include <atomic>
#include <thread>

#include "core/utilities/types.h"
#include "phare/phare.h"
#include "simulator/simulator.h"
#include "amr/wrappers/hierarchy.h"

#include "bench/core/bench.h"

#include "NumCpp.hpp"

namespace PHARE::real::bench::harris
{
std::size_t static constexpr dim    = 2;
std::size_t static constexpr interp = 1;
int constexpr cells                 = 100;
auto constexpr dl                   = .2;

using Param  = std::vector<double> const&;
using Return = std::shared_ptr<PHARE::core::Span<double>>;



template<typename T, typename SIZE = size_t>
class NCArraySpan : private core::StackVar<nc::NdArray<T>>, public core::Span<T, SIZE>
{
    using Vector = core::StackVar<nc::NdArray<T>>;
    using Span_  = core::Span<T, SIZE>;

public:
    using Vector::var;

    NCArraySpan(std::size_t size, T value = 0)
        : Vector{std::vector<T>(size, value)}
        , Span_{Vector::var.data(), Vector::var.size()}
    {
    }
};

using VecSpan = std::shared_ptr<NCArraySpan<double>>;


template<typename Vec>
auto asarray(Vec const& v0)
{
    return nc::asarray<double>(const_cast<std::vector<double>&>(v0), /*copy=*/false);
}

template<typename Ret = Return>
Ret bx(Param x, Param y)
{
    return std::make_shared<NCArraySpan<double>>(x.size(), 0);
}

template<typename Ret = Return>
Ret by(Param x, Param y)
{
    return std::make_shared<NCArraySpan<double>>(x.size(), 0);
}

template<typename Ret = Return>
Ret bz(Param x, Param y)
{
    return std::make_shared<NCArraySpan<double>>(x.size(), 0);
}

template<typename Ret = Return>
Ret b2(Param x, Param y)
{
    auto bx0  = bx<VecSpan>(x, y);
    auto& bx_ = bx0->var;

    auto by0  = by<VecSpan>(x, y);
    auto& by_ = bx0->var;

    auto bz0  = bz<VecSpan>(x, y);
    auto& bz_ = bx0->var;

    auto v0 = std::make_shared<NCArraySpan<double>>(x.size(), 0);
    auto& v = v0->var;

    assert(&v0->var[0] == &v[0]);
    v = nc::power(bx_, 2) + nc::power(by_, 2) + nc::power(bz_, 2);
    assert(&v0->var[0] == &v[0]);
    return v0;
}

template<typename Ret = Return>
Ret density(Param x0, Param y0)
{
    auto L  = cells;
    auto x  = asarray(x0);
    auto y  = asarray(y0);
    auto v0 = std::make_shared<core::VectorSpan<double>>(x0.size(), 0);
    auto v  = asarray(v0->var);
    assert(&v0->var[0] == &v[0]);

    KLOG(INF) << v0->var[10] << " " << v0->var.back();

    v = 0.2 + 1. / nc::power(nc::cosh((y - L * 0.3) / 0.5), 2)
        + 1. / nc::power(nc::cosh((y - L * 0.7) / 0.5), 2);

    assert(&v0->var[0] == &v[0]);

    KLOG(INF) << v0->var[10] << " " << v0->var.back();
    return v0;
}

Return T(Param x0, Param y0)
{
    auto rho0 = density<VecSpan>(x0, y0);
    auto rho  = asarray(rho0->var);
    auto b0   = b2<VecSpan>(x0, y0);
    auto b    = asarray(b0->var);

    auto v0 = std::make_shared<core::VectorSpan<double>>(x0.size(), 0);
    auto v  = asarray(v0->var);

    auto K0 = std::make_shared<core::VectorSpan<double>>(x0.size(), 1);
    auto K  = asarray(K0->var); // annoying

    assert(&v0->var[0] == &v[0]);
    v = 1. / rho * (K - b * 0.5);
    assert(&v0->var[0] == &v[0]);
    return v0;
}

Return vx(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x.size(), 0);
}

Return vy(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x.size(), 0);
}

Return vz(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x.size(), 0);
}

Return vthx(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vthy(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

Return vthz(Param x, Param y)
{
    return std::make_shared<core::VectorSpan<double>>(x);
}

PHARE::initializer::PHAREDict& createDict()
{
    using InitFunctionT = PHARE::initializer::InitFunction<dim>;

    auto& dict = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();

    auto& sim = dict["simulation"];

    sim["dimension"]            = 2;
    sim["interp_order"]         = 1;
    sim["refined_particle_nbr"] = 4;
    sim["time_step"]            = .001;
    sim["time_step_nbr"]        = 10;
    sim["boundary_type"]        = std::string{"periodic"};

    for (auto s : {"x", "y"})
    {
        sim["grid"]["nbr_cells"][s] = cells;
        sim["grid"]["meshsize"][s]  = dl;
        sim["grid"]["origin"][s]    = 0.;
    }

    auto& amr                              = sim["AMR"];
    amr["max_nbr_levels"]                  = 2;
    amr["nesting_buffer"]                  = std::vector<int>{0, 0};
    amr["refinement"]["tagging"]["method"] = std::string{"auto"};

    auto& algo                            = sim["algo"];
    algo["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};
    algo["ohm"]["resistivity"]            = .001;
    algo["ohm"]["hyper_resistivity"]      = .001;

    sim["ions"]["nbrPopulations"] = int{1};
    sim["ions"]["pop0"]["name"]   = std::string{"protons"};
    sim["ions"]["pop0"]["mass"]   = 1.;

    auto& pop_init                 = sim["ions"]["pop0"]["particle_initializer"];
    pop_init["name"]               = std::string{"maxwellian"};
    pop_init["nbr_part_per_cell"]  = int{100};
    pop_init["charge"]             = -1.;
    pop_init["basis"]              = std::string{"cartesian"};
    pop_init["density"]            = static_cast<InitFunctionT>(density<Return>);
    pop_init["bulk_velocity_x"]    = static_cast<InitFunctionT>(vx);
    pop_init["bulk_velocity_y"]    = static_cast<InitFunctionT>(vy);
    pop_init["bulk_velocity_z"]    = static_cast<InitFunctionT>(vz);
    pop_init["thermal_velocity_x"] = static_cast<InitFunctionT>(vthx);
    pop_init["thermal_velocity_y"] = static_cast<InitFunctionT>(vthy);
    pop_init["thermal_velocity_z"] = static_cast<InitFunctionT>(vthz);

    sim["electromag"]["name"]             = std::string{"EM"};
    sim["electromag"]["electric"]["name"] = std::string{"E"};
    sim["electromag"]["magnetic"]["name"] = std::string{"B"};

    auto& b_init          = sim["electromag"]["magnetic"]["initializer"];
    b_init["x_component"] = static_cast<InitFunctionT>(bx<Return>);
    b_init["y_component"] = static_cast<InitFunctionT>(by<Return>);
    b_init["z_component"] = static_cast<InitFunctionT>(bz<Return>);

    sim["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    sim["electrons"]["pressure_closure"]["Te"]   = 0.12;

    return dict;
}

} // namespace PHARE::real::bench::harris

#endif /*PHARE_BENCH_HARRIS_BENCH_H*/
