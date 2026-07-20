#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_USING_COMPUTED_FLUX_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_USING_COMPUTED_FLUX_HPP


#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_mhd_field_evolvers.hpp"

#include "core/inner_boundary/inner_bc_context.hpp"
#include "core/inner_boundary/inner_boundary_defs.hpp"

namespace PHARE::solver
{
template<typename MHDModel>
class EulerUsingComputedFlux
{
    using level_t                     = MHDModel::level_t;
    using gridlayout_type             = MHDModel::gridlayout_type;
    using state_type                  = MHDModel::state_type;
    using resources_manager_type      = MHDModel::resources_manager_type;
    using inner_boundary_manager_type = MHDModel::inner_boundary_manager_type;
    using Dispatchers_t               = Dispatchers<MHDModel>;

    using FiniteVolumeEuler_t = Dispatchers_t::FiniteVolumeEuler_t;
    using Faraday_t           = Dispatchers_t::Faraday_t;

public:
    EulerUsingComputedFlux() {}

    // we provide dt here because we sometimes need it to be different from newTime-currentTime, for
    // example in the case of some rk integration methods
    void operator()(MHDModel& model, auto& state, auto& statenew, auto& E, auto& fluxes, auto& bc,
                    level_t& level, double const newTime, double const dt)
    {
        FiniteVolumeEuler_t{level, model}(newTime, state, statenew, fluxes, dt);

        Faraday_t{level, model}(state.B, E, statenew.B, dt);

        resources_manager_type& rm = *model.resourcesManager;

        // Pin the inner-boundary inactive region to a canonical safe physical state.
        // We cannot simply copy the previous-step values into `statenew` because some
        // integrators (e.g. TVDRK3) reuse a single buffer for both `state` and `statenew`;
        // in that aliased case the copy is a self-assignment and the drift from the euler
        // update / Faraday persists. Recomputing the safe state every substep is both
        // alias-safe and matches what `mhd_level_initializer` does at init.
        if (model.hasInnerBoundary())
        {
            inner_boundary_manager_type& ibm = *model.innerBoundaryManager;

            amr::visitLevel<gridlayout_type>(
                level, rm,
                [&](auto& layout, auto&&, auto&&) {
                    // Faces (B) first, then the cell-centered conserved moments.
                    ibm.setSafeState(statenew.B, layout);
                    ibm.setSafeState(statenew.rho, layout);
                    ibm.setSafeState(statenew.rhoV, layout);
                    ibm.setSafeState(statenew.Etot, layout);
                },
                ibm, statenew);
        }

        bc.fillMagneticGhosts(statenew.B, level, newTime);

        bc.fillMomentsGhosts(statenew, level, newTime);

        // Application of inner boundary conditions on the moments (priority-ordered).
        if (model.hasInnerBoundary())
        {
            inner_boundary_manager_type& ibm = *model.innerBoundaryManager;
            core::InnerBCContext<state_type> ctx{statenew, state, newTime, dt};

            amr::visitLevel<gridlayout_type>(
                level, rm,
                [&](auto& layout, auto&&, auto&&) { ibm.applyToMoments(layout, ctx); }, ibm, state,
                statenew);
        }
    }
};


} // namespace PHARE::solver

#endif
