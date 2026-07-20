#ifndef PHARE_SOLVER_TIME_INTEGRATOR_HPP
#define PHARE_SOLVER_TIME_INTEGRATOR_HPP

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "initializer/data_provider.hpp"
#include "phare_simulator_options.hpp"

#include "amr/solvers/time_integrator/euler_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk2_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk3_integrator.hpp"
#include "amr/solvers/time_integrator/ssprk4_5_integrator.hpp"

namespace PHARE::solver
{

template<typename FVMethodStrategy, typename MHDModel,
         typename MessengerT = amr::MHDMessenger<MHDModel>>
class TimeIntegrator
{
    using Base_t = BaseMHDTimestepper<MHDModel, MessengerT>;

    static std::unique_ptr<Base_t> make_(MHDOpts::TimeIntegratorType t,
                                         PHARE::initializer::PHAREDict const& dict)
    {
        switch (t)
        {
            case MHDOpts::TimeIntegratorType::Euler:
                return std::make_unique<EulerIntegrator<FVMethodStrategy, MHDModel, MessengerT>>(
                    dict);
            case MHDOpts::TimeIntegratorType::TVDRK2:
                return std::make_unique<TVDRK2Integrator<FVMethodStrategy, MHDModel, MessengerT>>(
                    dict);
            case MHDOpts::TimeIntegratorType::TVDRK3:
                return std::make_unique<TVDRK3Integrator<FVMethodStrategy, MHDModel, MessengerT>>(
                    dict);
            case MHDOpts::TimeIntegratorType::SSPRK4_5:
                return std::make_unique<SSPRK4_5Integrator<FVMethodStrategy, MHDModel, MessengerT>>(
                    dict);
            default: throw std::runtime_error("TimeIntegrator: unsupported time integrator type");
        }
    }

    std::unique_ptr<Base_t> impl_;

public:
    using Messenger = MessengerT;

    TimeIntegrator(PHARE::initializer::PHAREDict const& dict)
        : impl_{make_(core::fromString<MHDOpts::TimeIntegratorType>(
                          cppdict::get_value(dict, "time_integrator_type", std::string{"TVDRK3"})),
                      dict)}
    {
    }

    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& bc, MHDModel::level_t& level,
                    double const currentTime, double const newTime)
    {
        (*impl_)(model, state, fluxes, bc, level, currentTime, newTime);
    }

    void registerResources(MHDModel& model) { impl_->registerResources(model); }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        impl_->allocate(model, patch, allocateTime);
    }

    void fillMessengerInfo(auto& info) const { impl_->fillMessengerInfo(info); }

    auto exposeFluxes() { return impl_->exposeFluxes(); }

    auto exposeFluxes() const { return std::as_const(*impl_).exposeFluxes(); }
};

} // namespace PHARE::solver

#endif // PHARE_SOLVER_TIME_INTEGRATOR_HPP
