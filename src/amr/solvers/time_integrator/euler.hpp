#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/compute_fluxes_and_sources.hpp"
#include "amr/solvers/time_integrator/euler_using_computed_flux.hpp"

namespace PHARE::solver
{
template<typename FVMethodStrategy, typename MHDModel>
class Euler
{
    using level_t                  = MHDModel::level_t;
    using ComputeFluxesAndSources_t = ComputeFluxesAndSources<FVMethodStrategy, MHDModel>;
    using EulerUsingComputedFlux_t = EulerUsingComputedFlux<MHDModel>;

public:
    Euler(PHARE::initializer::PHAREDict const& dict)
        : compute_fluxes_and_sources_{dict}
    {
    }

    void operator()(MHDModel& model, auto& state, auto& statenew, auto& fluxes, auto& sources,
                    auto& bc, level_t& level, double const currentTime, double const newTime,
                    double dt = std::nan(""))
    {
        if (std::isnan(dt))
            dt = newTime - currentTime;

        compute_fluxes_and_sources_(model, state, fluxes, sources, bc, level, newTime);

        euler_using_computed_flux_(model, state, statenew, state.E, fluxes, sources, bc, level,
                                   newTime, dt);
    }

    void registerResources(MHDModel& model) { compute_fluxes_and_sources_.registerResources(model); }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        compute_fluxes_and_sources_.allocate(model, patch, allocateTime);
    }

private:
    ComputeFluxesAndSources_t compute_fluxes_and_sources_;
    EulerUsingComputedFlux_t euler_using_computed_flux_;
};
} // namespace PHARE::solver

#endif
