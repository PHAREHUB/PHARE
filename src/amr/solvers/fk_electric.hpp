#ifndef PHARE_AMR_SOLVERS_FK_ELECTRIC_HPP
#define PHARE_AMR_SOLVERS_FK_ELECTRIC_HPP

#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/numerics/ampere/ampere_maxwell.hpp"
#include "core/numerics/esirkepov/esirkepov.hpp"
#include "core/numerics/filter/binomial_filter.hpp"
#include "core/numerics/marder/marder.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "mpi/mpi_utils.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"

#include <SAMRAI/hier/HierarchyNeighbors.h>
#include <SAMRAI/hier/PatchHierarchy.h>

#include <string>
#include <unordered_map>
#include <vector>


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



/**
 * @brief regionEdgeNodes returns the primal nodes of a patch (AMR indexes of its `cells`, node i
 * the lower corner of cell i) whose 2^dim cells are some in `region` and some not, one box per
 * node. Such a node is on a face of a region box: only those faces, within the patch, are
 * searched, and a node on several faces is listed once per face.
 *
 * With boundary_cleaning_mask the cleaning skips these nodes (Marder, F = 0 there):
 * - region = the patch and its same-level neighbours (nonLevelGhostBox): the edge of a refined
 *   level, whose residual reads the coarse E ghost (t5_fine_gauss.md, F0; pV);
 * - region = the cells the finer level covers: the edge of the synchronized E on a coarser level
 *   (t10_interface_model.py, coarse_edge_mask).
 * Inner corners of a level are on its edge too, though their residual only reads level E: a
 * superset, as for the finer boxes, whatever edges the synchronization writes on their faces.
 */
template<typename Box_t>
auto regionEdgeNodes(Box_t const& cells, std::vector<Box_t> const& region)
{
    constexpr auto dimension = Box_t::dimension;
    constexpr auto corners   = std::size_t{1} << dimension;

    auto const nodeBox = [](Box_t box) {
        for (std::size_t d = 0; d < dimension; ++d)
            box.upper[d] += 1;
        return box;
    };
    auto const patchNodes = nodeBox(cells);

    // the 2^dim cells around node i are i - 1 and i in each direction
    auto const onEdge = [&](auto const& node) {
        std::size_t in = 0;
        for (std::size_t corner = 0; corner < corners; ++corner)
        {
            auto cell = node;
            for (std::size_t d = 0; d < dimension; ++d)
                cell[d] -= static_cast<int>((corner >> d) & 1u);
            in += core::isIn(cell, region);
        }
        return in > 0 and in < corners;
    };

    std::vector<Box_t> nodes;
    for (auto const& box : region)
        for (std::size_t d = 0; d < dimension; ++d)
            for (auto const face : {box.lower[d], box.upper[d] + 1})
            {
                auto faceNodes     = nodeBox(box);
                faceNodes.lower[d] = faceNodes.upper[d] = face;
                if (auto const searched = faceNodes * patchNodes)
                    core::for_each_point(*searched, [&](auto const& node) {
                        if (onEdge(node))
                            nodes.emplace_back(node, node);
                    });
            }
    return nodes;
}



template<std::size_t dimension>
using PatchBoxes = std::unordered_map<std::string /*patch id*/, std::vector<core::Box<int, dimension>>>;

/**
 * @brief finerLevelCoverage returns, for every local patch of level ln, the boxes of the next
 * finer level near it (SAMRAI's coarse to fine connector, width 1: boxes of other ranks and
 * periodic images included), coarsened to level ln: the cells the finer level synchronizes.
 * Empty without a finer level. Collective (the connector may be built).
 */
template<std::size_t dimension>
auto finerLevelCoverage(SAMRAI::hier::PatchHierarchy const& hierarchy, int const ln)
{
    PatchBoxes<dimension> covered;
    if (!hierarchy.finerLevelExists(ln))
        return covered;

    auto const ratio = hierarchy.getPatchLevel(ln + 1)->getRatioToCoarserLevel();
    SAMRAI::hier::HierarchyNeighbors const neighbors{hierarchy, ln, ln + 1,
                                                     /*do_same_level_nbrs=*/false};
    for (auto const& patch : *hierarchy.getPatchLevel(ln))
    {
        auto& boxes = covered[amr::to_string(patch->getGlobalId())];
        for (auto fine : neighbors.getFinerLevelNeighbors(patch->getBox(), ln))
        {
            fine.coarsen(ratio);
            boxes.push_back(amr::phare_box_from<dimension>(fine));
        }
    }
    return covered;
}

/**
 * @brief cleaningSkippedNodes returns, for every local patch of level ln with any, the nodes the
 * cleaning skips (regionEdgeNodes): on a refined level the edge of the level, and the edge of the
 * `covered` cells of the patch (finerLevelCoverage, of the last synchronization).
 */
template<typename GridLayout>
auto cleaningSkippedNodes(SAMRAI::hier::PatchHierarchy const& hierarchy, int const ln,
                          PatchBoxes<GridLayout::dimension> const& covered)
{
    constexpr auto dimension = GridLayout::dimension;
    auto const& level        = *hierarchy.getPatchLevel(ln);

    PatchBoxes<dimension> skipped;
    auto const add = [&](auto const& patch, auto const& region) {
        auto nodes = regionEdgeNodes(amr::phare_box_from<dimension>(patch.getBox()), region);
        if (!nodes.empty())
        {
            auto& all = skipped[amr::to_string(patch.getGlobalId())];
            all.insert(all.end(), nodes.begin(), nodes.end());
        }
    };

    if (ln > 0)
    {
        auto const neighbors = amr::makeSameLevelNeighbors(hierarchy, ln);
        for (auto const& patch : level)
            add(*patch, amr::makeNonLevelGhostBoxFor<GridLayout>(*patch, neighbors));
    }
    for (auto const& patch : level)
        if (auto const it = covered.find(amr::to_string(patch->getGlobalId())); it != covered.end())
            add(*patch, it->second);
    return skipped;
}



/**
 * @brief advanceFKElectric is steps 1 and 2 of SolverFK: J and rho from the particles at time
 * (depositFKCurrent), then E, at time - dt/2, cleaned (Marder, d from model.divergence_cleaning on
 * this level's mesh) and advanced to time + dt/2 (Ampere-Maxwell). E is the model E, or the copy of
 * it the solver advances one step ahead for a finer level (P4 plan D1). E ghosts are read, not
 * filled; labels are the caller's. The cleaning of a patch skips its nodes in `skipped`
 * (cleaningSkippedNodes; empty without boundary_cleaning_mask).
 */
template<typename FKModel, typename Messenger, typename VecField>
void advanceFKElectric(typename FKModel::level_t& level, FKModel& model, Messenger& messenger,
                       VecField& E, double const dt, double const time,
                       PatchBoxes<FKModel::dimension> const& skipped)
{
    using GridLayout = FKModel::gridlayout_type;

    auto& rm  = *model.resourcesManager;
    auto& B   = model.state.electromag.B;
    auto& J   = model.state.J;
    auto& rho = model.state.rho;

    depositFKCurrent(level, model, messenger, dt, /*back=*/1., /*front=*/1., time);

    if (model.divergence_cleaning > 0.)
        for (auto& patch : rm.enumerate(level, rho, E))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            auto const d      = core::marderCoefficient(layout, model.divergence_cleaning, dt);
            core::Marder const marder{layout, model.c2(), d};
            auto const it = skipped.empty() ? skipped.end()
                                            : skipped.find(amr::to_string(patch->getGlobalId()));
            if (it == skipped.end())
                marder(E, rho, dt);
            else
                marder(E, rho, dt, it->second);
        }
    for (auto& patch : rm.enumerate(level, B, J, E))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        core::AmpereMaxwell{layout, model.c2()}(B, J, E, dt);
    }
}


/**
 * @brief cleanElectric applies `passes` Marder steps of strength s to E on a level whose rho
 * (depositFKCurrent) stays fixed: the E of a level (re)created during the run, before its first
 * step (P5a). The same-level E ghosts are filled again before every pass; the level ghosts keep
 * what the caller filled, the coarser E interpolated at the time of E, which a pass does not
 * change. Each pass multiplies the Gauss residual mode k by
 * 1 - s sum_i (4/dl_i^2) sin^2(k_i dl_i / 2) / sum_i (2/dl_i^2) away from the level boundary.
 * The passes on a patch skip its nodes in `skipped` (cleaningSkippedNodes; empty without
 * boundary_cleaning_mask).
 */
template<typename FKModel, typename Messenger, typename VecField>
void cleanElectric(typename FKModel::level_t& level, FKModel& model, Messenger& messenger,
                   VecField& E, double const dt, double const s, std::size_t const passes,
                   double const time, PatchBoxes<FKModel::dimension> const& skipped)
{
    using GridLayout = FKModel::gridlayout_type;

    auto& rm  = *model.resourcesManager;
    auto& rho = model.state.rho;
    for (std::size_t pass = 0; pass < passes; ++pass)
    {
        messenger.fillElectricPatchGhosts(E, level, time);
        for (auto& patch : rm.enumerate(level, rho, E))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            auto const d      = core::marderCoefficient(layout, s, dt);
            core::Marder const marder{layout, model.c2(), d};
            auto const it = skipped.empty() ? skipped.end()
                                            : skipped.find(amr::to_string(patch->getGlobalId()));
            if (it == skipped.end())
                marder(E, rho, dt);
            else
                marder(E, rho, dt, it->second);
        }
    }
}

} // namespace PHARE::solver

#endif
