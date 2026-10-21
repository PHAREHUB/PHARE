#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_USING_COMPUTED_FLUX_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_USING_COMPUTED_FLUX_HPP


#include "amr/solvers/solver_mhd_field_evolvers.hpp"

namespace PHARE::solver
{
template<typename MHDModel>
class EulerUsingComputedFlux
{
    using level_t = MHDModel::level_t;
    // using Layout        = MHDModel::gridlayout_type;
    using Dispatchers_t = Dispatchers<MHDModel>;

    using FiniteVolumeEuler_t = Dispatchers_t::FiniteVolumeEuler_t;
    using Faraday_t           = Dispatchers_t::Faraday_t;
    using GridLayout          = MHDModel::gridlayout_type;

public:
    EulerUsingComputedFlux() {}

    // we provide dt here because we sometimes need it to be different from newTime-currentTime, for
    // example in the case of some rk integration methods
    void operator()(MHDModel& model, auto& state, auto& statenew, auto& E, auto& fluxes,
                    auto& sources, auto& bc, level_t& level, double const newTime, double const dt)
    {
        FiniteVolumeEuler_t{level, model}(state, statenew, fluxes, dt);
        TimeSetter{level, model, newTime}(state.rho, state.rhoV, state.Etot1);

        Faraday_t{level, model}(state.B1, E, statenew.B1, dt);

        applySources_(model, statenew, sources, level, dt);

        TimeSetter{level, model, newTime}(statenew.B1);

        bc.fillMagneticGhosts(statenew.B1, level, newTime);

        bc.fillMomentsGhosts(statenew, level, newTime);
    }

private:
    void applySources_(MHDModel& model, auto& statenew, auto& sources, level_t& level,
                       double const dt)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, statenew, sources))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);

            auto& Etot1    = statenew.Etot1;
            auto const& sE = sources.Etot_source;
            layout.evalOnBox(Etot1, [&](auto&... args) { Etot1(args...) += dt * sE(args...); });

            for (auto const c : {core::Component::X, core::Component::Y, core::Component::Z})
            {
                auto& B1c       = statenew.B1(c);
                auto const& sBc = sources.B1_source(c);
                layout.evalOnBox(B1c, [&](auto&... args) { B1c(args...) += dt * sBc(args...); });
            }
        }
    }
};


} // namespace PHARE::solver

#endif
