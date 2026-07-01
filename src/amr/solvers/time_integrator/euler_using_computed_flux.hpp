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
        FiniteVolumeEuler_t{level, model}(newTime, state, statenew, fluxes, dt);

        // Faraday evolves only the perturbation B1; the background B0 lives in the model and is
        // shared by every stage (no per-stage copy).
        Faraday_t{level, model}(state.B1, E, statenew.B1, dt);

        // Apply the body sources (dU/dt = -div F + S): statenew += dt * S. `sources` holds either a
        // single stage's sources or the RK-accumulated Butcher sources, consistently with `fluxes`
        // and `E`. Zero for a static background B0, so this is an exact no-op then.
        applySources_(model, statenew, sources, level, dt);

        bc.fillMagneticGhosts(statenew.B1, level, newTime);

        bc.fillMomentsGhosts(statenew, level, newTime);
    }

private:
    // statenew.Etot1 += dt * Etot_source (cell-centered) and statenew.B1 += dt * B1_source (face-
    // centered), on the same domain the finite-volume / Faraday updates just wrote.
    void applySources_(MHDModel& model, auto& statenew, auto& sources, level_t& level,
                       double const dt)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, statenew, sources))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);

            auto& Etot1       = statenew.Etot1;
            auto const& es    = sources.Etot_source;
            layout.evalOnBox(Etot1,
                             [&](auto&... args) mutable { Etot1(args...) += dt * es(args...); });

            for (auto const& c : {core::Component::X, core::Component::Y, core::Component::Z})
            {
                auto& B1c      = statenew.B1(c);
                auto const& bs = sources.B1_source(c);
                layout.evalOnBox(B1c,
                                 [&](auto&... args) mutable { B1c(args...) += dt * bs(args...); });
            }
        }
    }
};


} // namespace PHARE::solver

#endif
