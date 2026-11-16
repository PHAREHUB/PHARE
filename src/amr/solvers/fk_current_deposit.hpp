#ifndef PHARE_AMR_SOLVERS_FK_CURRENT_DEPOSIT_HPP
#define PHARE_AMR_SOLVERS_FK_CURRENT_DEPOSIT_HPP

#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/numerics/esirkepov/esirkepov.hpp"
#include "mpi/mpi_utils.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"


namespace PHARE::solver
{
/**
 * @brief depositFKCurrent sets the fully kinetic model current J on a level: the Esirkepov
 * current of every owned domain particle, over its segment x - back (dt/2) v -> x + front (dt/2) v,
 * summed on patch borders. Domain particles are the owned ones: the push (or the level
 * initializer) erased the others. Patch ghost particles never deposit, their owner does.
 * A particle crossing a cell in a step stops the run on all ranks, the error names its population.
 */
template<typename FKModel, typename Messenger>
void depositFKCurrent(typename FKModel::level_t& level, FKModel& model, Messenger& messenger,
                      double const dt, double const back, double const front, double const time)
{
    PHARE_LOG_SCOPE(1, "depositFKCurrent");

    using GridLayout = FKModel::gridlayout_type;
    core::Esirkepov<FKModel::dimension, FKModel::interp_order> esirkepov;

    auto& rm   = *model.resourcesManager;
    auto& ions = model.state.ions;
    auto& J    = model.state.J;

    try
    {
        for (auto& patch : rm.enumerate(level, ions, J))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            J.zero();
            for (auto const& pop : ions)
            {
                try
                {
                    esirkepov(pop.domainParticles(), J, layout, pop.charge(), dt, back, front);
                }
                catch (core::DictionaryException& ex)
                {
                    throw ex("population", pop.name());
                }
            }
        }
    }
    catch (core::DictionaryException const& ex)
    {
        PHARE_LOG_ERROR(ex());
    }
    if (mpi::any_errors())
        throw core::DictionaryException{}("ID", "depositFKCurrent");

    TimeSetter{level, model, time}(J);
    messenger.fillCurrentBorders(J, level, time);
}

} // namespace PHARE::solver

#endif
