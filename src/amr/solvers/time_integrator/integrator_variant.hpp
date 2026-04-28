#ifndef PHARE_CORE_NUMERICS_INTEGRATOR_VARIANT_HPP
#define PHARE_CORE_NUMERICS_INTEGRATOR_VARIANT_HPP

#include <stdexcept>
#include <variant>

#include "initializer/data_provider.hpp"
#include "phare_simulator_options.hpp"

#include "amr/solvers/time_integrator/euler_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk2_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk3_integrator.hpp"
#include "amr/solvers/time_integrator/ssprk4_5_integrator.hpp"

namespace PHARE::solver
{

namespace detail
{
inline MHDOpts::TimeIntegratorType parse_time_integrator_type(std::string const& s)
{
    if (s == "Euler" || s == "euler") return MHDOpts::TimeIntegratorType::Euler;
    if (s == "TVDRK2" || s == "tvdrk2") return MHDOpts::TimeIntegratorType::TVDRK2;
    if (s == "TVDRK3" || s == "tvdrk3") return MHDOpts::TimeIntegratorType::TVDRK3;
    if (s == "SSPRK4_5" || s == "ssprk4_5") return MHDOpts::TimeIntegratorType::SSPRK4_5;
    throw std::runtime_error("Unknown time integrator type: " + s);
}
} // namespace detail

template<typename FVMethodStrategy, typename MHDModel>
class IntegratorVariant
{
    using Euler_t    = EulerIntegrator<FVMethodStrategy, MHDModel>;
    using TVDRK2_t   = TVDRK2Integrator<FVMethodStrategy, MHDModel>;
    using TVDRK3_t   = TVDRK3Integrator<FVMethodStrategy, MHDModel>;
    using SSPRK4_5_t = SSPRK4_5Integrator<FVMethodStrategy, MHDModel>;

    using Variant_t = std::variant<Euler_t, TVDRK2_t, TVDRK3_t, SSPRK4_5_t>;

    static Variant_t make_(MHDOpts::TimeIntegratorType t,
                           PHARE::initializer::PHAREDict const& dict)
    {
        switch (t)
        {
            case MHDOpts::TimeIntegratorType::Euler:   return Euler_t{dict};
            case MHDOpts::TimeIntegratorType::TVDRK2:  return TVDRK2_t{dict};
            case MHDOpts::TimeIntegratorType::TVDRK3:  return TVDRK3_t{dict};
            case MHDOpts::TimeIntegratorType::SSPRK4_5: return SSPRK4_5_t{dict};
            default:
                throw std::runtime_error("IntegratorVariant: unsupported time integrator type");
        }
    }

    Variant_t v_;

public:
    IntegratorVariant(PHARE::initializer::PHAREDict const& dict)
        : v_{make_(
              detail::parse_time_integrator_type(
                  cppdict::get_value(dict, "time_integrator_type", std::string{"TVDRK3"})),
              dict)}
    {
    }

    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& bc,
                    typename MHDModel::level_t& level, double const currentTime,
                    double const newTime)
    {
        std::visit(
            [&](auto& e) { e(model, state, fluxes, bc, level, currentTime, newTime); }, v_);
    }

    void registerResources(MHDModel& model)
    {
        std::visit([&](auto& e) { e.registerResources(model); }, v_);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        std::visit([&](auto const& e) { e.allocate(model, patch, allocateTime); }, v_);
    }

    void fillMessengerInfo(auto& info) const
    {
        std::visit([&](auto const& e) { e.fillMessengerInfo(info); }, v_);
    }

    auto exposeFluxes()
    {
        return std::visit([](auto& e) { return e.exposeFluxes(); }, v_);
    }

    auto exposeFluxes() const
    {
        return std::visit([](auto const& e) { return e.exposeFluxes(); }, v_);
    }
};

} // namespace PHARE::solver

#endif // PHARE_CORE_NUMERICS_INTEGRATOR_VARIANT_HPP
