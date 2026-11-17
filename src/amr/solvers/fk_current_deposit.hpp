#ifndef PHARE_AMR_SOLVERS_FK_CURRENT_DEPOSIT_HPP
#define PHARE_AMR_SOLVERS_FK_CURRENT_DEPOSIT_HPP

#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/numerics/esirkepov/esirkepov.hpp"
#include "core/numerics/filter/binomial_filter.hpp"
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
 *
 * With model.current_filter_passes > 0, J is then filtered by binomial passes. A pass reads one
 * ghost layer and writes the physical nodes only, so after a pass the ghosts still hold the
 * unfiltered neighbor values: the ghosts are filled again before every pass, fill then filter,
 * never two passes on one fill. (A ghost width of 2, order 1, could serve 2 passes per fill:
 * a later optimization keeping this loop.) The solver and the level initializer both deposit
 * here, so the E half step back of t=0 sees the filtered current too.
 *
 * With model.divergence_cleaning > 0, the same particles also deposit the charge density rho at
 * the start of their segment, x - back (dt/2) v: the charge Gauss's law pairs with the E before
 * this J updates it (E^{n-1/2} in the solver, E^{-1/2} in the initializer). rho is summed on
 * patch borders and filtered like J, the filter commuting with the Yee divergence. Without
 * cleaning rho is not deposited.
 *
 * On a refined level, the level ghost particles (the pushable copies each patch holds in its
 * level ghost band) then add their J and rho, after the border sums: a patch deposits its own
 * copies on its own nodes, physical and ghost, so a node shared by several patches gets each
 * level ghost particle once, and the band nodes are complete. The ghost fills of the filter
 * copy from same-level patches only, the band keeps that deposit. With 2 passes or more, the
 * band nodes are not filtered again after the first pass: S^N is not exact on the last N-1
 * node layers next to the level boundary, a Gauss residual source the cleaning absorbs.
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
    auto& rho  = model.state.rho;

    bool const cleaning = model.divergence_cleaning > 0.;

    // J (and rho) of the given particles of every population, added to what the patch holds
    auto const deposit = [&](auto&& particlesOf) {
        try
        {
            for (auto& patch : rm.enumerate(level, ions, J, rho))
            {
                auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
                for (auto const& pop : ions)
                {
                    auto const& particles = particlesOf(pop);
                    try
                    {
                        esirkepov(particles, J, layout, pop.charge(), dt, back, front);
                    }
                    catch (core::DictionaryException& ex)
                    {
                        throw ex("population", pop.name());
                    }
                    if (cleaning)
                        esirkepov.depositCharge(particles, rho, layout, pop.charge(),
                                                -back * 0.5 * dt);
                }
            }
        }
        catch (core::DictionaryException const& ex)
        {
            PHARE_LOG_ERROR(ex());
        }
        if (mpi::any_errors())
            throw core::DictionaryException{}("ID", "depositFKCurrent");
    };

    for (auto& patch : rm.enumerate(level, J, rho))
    {
        J.zero();
        if (cleaning)
            rho.zero();
    }
    deposit([](auto const& pop) -> auto const& { return pop.domainParticles(); });

    // J at the time of the particles, rho at its own, the start of their segments
    auto const rhoTime = time - back * 0.5 * dt;
    TimeSetter{level, model, time}(J);
    messenger.fillCurrentBorders(J, level, time);
    if (cleaning)
    {
        TimeSetter{level, model, rhoTime}(rho);
        messenger.fillChargeDensityBorders(rho, level, rhoTime);
    }

    // after the border sums: every patch holding a copy of a level ghost particle deposits it
    // on its own nodes, so a node shared by patches gets it once from each (P4 plan D5)
    if (level.getLevelNumber() > 0)
        deposit([](auto const& pop) -> auto const& { return pop.levelGhostParticles(); });

    // fill then filter, every pass
    auto const filter = [&](auto& field, auto&& fillGhosts) {
        for (std::size_t pass = 0; pass < model.current_filter_passes; ++pass)
        {
            fillGhosts();
            for (auto& patch : rm.enumerate(level, field))
                core::BinomialFilter{amr::layoutFromPatch<GridLayout>(*patch)}(field);
        }
    };

    filter(J, [&]() { messenger.fillCurrentGhosts(J, level, time); });
    if (cleaning)
        filter(rho, [&]() { messenger.fillChargeDensityGhosts(rho, level, rhoTime); });
}

} // namespace PHARE::solver

#endif
