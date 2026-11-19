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

#include <SAMRAI/hier/BoxLevel.h>
#include <SAMRAI/hier/HierarchyNeighbors.h>
#include <SAMRAI/hier/PatchHierarchy.h>

#include <array>
#include <cstdint>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>


namespace PHARE::solver
{
// deposit(particles, layout, charge) for the particles particlesOf(pop) of every population on
// every patch of the level; an error of a deposit stops the run on all ranks, naming its population
template<typename FKModel>
void depositOnLevel(typename FKModel::level_t& level, FKModel& model, auto&& particlesOf,
                    auto&& deposit)
{
    using GridLayout = FKModel::gridlayout_type;
    auto& ions       = model.state.ions;
    try
    {
        for (auto& patch :
             model.resourcesManager->enumerate(level, ions, model.state.J, model.state.rho))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            for (auto const& pop : ions)
                try
                {
                    deposit(particlesOf(pop), layout, pop.charge());
                }
                catch (core::DictionaryException& ex)
                {
                    throw ex("population", pop.name());
                }
        }
    }
    catch (core::DictionaryException const& ex)
    {
        PHARE_LOG_ERROR(ex());
    }
    if (mpi::any_errors())
        throw core::DictionaryException{}("ID", "depositFKCurrent");
}


// fill then filter `field` on the level, model.current_filter_passes times
template<typename FKModel>
void filterOnLevel(typename FKModel::level_t& level, FKModel& model, auto& field, auto&& fillGhosts)
{
    using GridLayout = FKModel::gridlayout_type;
    for (std::size_t pass = 0; pass < model.current_filter_passes; ++pass)
    {
        fillGhosts();
        for (auto& patch : model.resourcesManager->enumerate(level, field))
            core::BinomialFilter{amr::layoutFromPatch<GridLayout>(*patch)}(field);
    }
}


/**
 * @brief depositFKCharge sets the fully kinetic model charge density rho on a level at time: that
 * of every owned domain particle at x + shift (dt/2) v, summed on patch borders, then on a refined
 * level that of the level ghost particles (as their current, depositFKCurrent), filtered like J,
 * the filter commuting with the Yee divergence. On a refined level its ghost nodes then take their
 * owners' values: the filter leaves them unfiltered, and the band nodes of a patch hold its own
 * band only. Gauss's law on the edge nodes of the level reads them (gaussOnLevelBorder).
 */
template<typename FKModel, typename Messenger>
void depositFKCharge(typename FKModel::level_t& level, FKModel& model, Messenger& messenger,
                     double const dt, double const shift, double const time)
{
    PHARE_LOG_SCOPE(1, "depositFKCharge");

    core::Esirkepov<FKModel::dimension, FKModel::interp_order> esirkepov;
    auto& rho = model.state.rho;

    auto const deposit = [&](auto const& particles, auto const& layout, double const charge) {
        esirkepov.depositCharge(particles, rho, layout, charge, shift * 0.5 * dt);
    };

    for (auto& patch : model.resourcesManager->enumerate(level, rho))
        rho.zero();
    depositOnLevel(
        level, model, [](auto const& pop) -> auto const& { return pop.domainParticles(); },
        deposit);
    TimeSetter{level, model, time}(rho);
    messenger.fillChargeDensityBorders(rho, level, time);

    if (level.getLevelNumber() > 0)
    {
        depositOnLevel(
            level, model, [](auto const& pop) -> auto const& { return pop.levelGhostParticles(); },
            deposit);
        messenger.equalizeChargeDensityBorders(rho, level, time);
    }

    filterOnLevel(level, model, rho,
                  [&]() { messenger.fillChargeDensityGhosts(rho, level, time); });
    if (level.getLevelNumber() > 0)
        messenger.fillChargeDensityGhosts(rho, level, time);
}


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
 * With model.divergence_cleaning > 0, and on a refined level, the same particles also deposit the
 * charge density rho at the start of their segment, x - back (dt/2) v (depositFKCharge): the
 * charge Gauss's law pairs with the E before this J updates it (E^{n-1/2} in the solver, E^{-1/2}
 * in the initializer). The cleaning reads it, and on a refined level the faces of E outside it
 * (gaussOnLevelBorder). Otherwise rho is not deposited.
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

    core::Esirkepov<FKModel::dimension, FKModel::interp_order> esirkepov;
    auto& J = model.state.J;

    auto const deposit = [&](auto const& particles, auto const& layout, double const charge) {
        esirkepov(particles, J, layout, charge, dt, back, front);
    };

    for (auto& patch : model.resourcesManager->enumerate(level, J))
        J.zero();
    depositOnLevel(
        level, model, [](auto const& pop) -> auto const& { return pop.domainParticles(); },
        deposit);
    TimeSetter{level, model, time}(J);
    messenger.fillCurrentBorders(J, level, time);

    // after the border sums: every patch holding a copy of a level ghost particle deposits it
    // on its own nodes, so a node shared by patches gets it once from each (P4 plan D5). Each
    // patch adds these terms in its own order: the max makes the copies equal again
    if (level.getLevelNumber() > 0)
    {
        depositOnLevel(
            level, model, [](auto const& pop) -> auto const& { return pop.levelGhostParticles(); },
            deposit);
        messenger.equalizeCurrentBorders(J, level, time);
    }

    filterOnLevel(level, model, J, [&]() { messenger.fillCurrentGhosts(J, level, time); });

    if (model.divergence_cleaning > 0. or level.getLevelNumber() > 0)
        depositFKCharge(level, model, messenger, dt, -back, time - back * 0.5 * dt);
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
 * @brief finerLevelCoverage returns the cells of level ln the next finer level covers, the cells
 * it synchronizes: every box of the finer level, on every rank (SAMRAI's globalized box level),
 * coarsened to level ln, and its images one period away in each periodic direction. It does not
 * depend on how level ln is cut into patches: it outlives a new level ln object (a rebalanced
 * root, P5b). Empty without a finer level. Collective (the globalized box level may be built).
 */
template<std::size_t dimension>
auto finerLevelCoverage(SAMRAI::hier::PatchHierarchy const& hierarchy, int const ln)
{
    std::vector<core::Box<int, dimension>> covered;
    if (!hierarchy.finerLevelExists(ln))
        return covered;

    auto const& finer = *hierarchy.getPatchLevel(ln + 1);
    auto const ratio  = finer.getRatioToCoarserLevel();
    auto const period = hierarchy.getGridGeometry()->getPeriodicShift(
        hierarchy.getPatchLevel(ln)->getRatioToLevelZero());

    // no shift, and -1 or +1 period along every periodic direction (0 along the others)
    std::vector<std::array<int, dimension>> shifts(1);
    for (std::size_t d = 0; d < dimension; ++d)
        if (period[d] != 0)
        {
            auto const unshifted = shifts;
            for (auto const& base : unshifted)
                for (int const k : {-1, 1})
                {
                    auto shift = base;
                    shift[d]   = k * period[d];
                    shifts.push_back(shift);
                }
        }

    for (auto fine : finer.getBoxLevel()->getGlobalizedVersion().getGlobalBoxes())
    {
        if (fine.isPeriodicImage())
            continue;
        fine.coarsen(ratio);
        auto const box = amr::phare_box_from<dimension>(fine);
        for (auto const& shift : shifts)
        {
            auto image = box;
            for (std::size_t d = 0; d < dimension; ++d)
            {
                image.lower[d] += shift[d];
                image.upper[d] += shift[d];
            }
            covered.push_back(image);
        }
    }
    return covered;
}

/**
 * @brief cleaningSkippedNodes returns, for every local patch of level ln with any, the nodes the
 * cleaning skips (regionEdgeNodes): on a refined level the edge of the level, and in the patch the
 * edge of the `covered` cells (finerLevelCoverage, of the last synchronization).
 */
template<typename GridLayout>
auto cleaningSkippedNodes(SAMRAI::hier::PatchHierarchy const& hierarchy, int const ln,
                          std::vector<core::Box<int, GridLayout::dimension>> const& covered)
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
    // only the covered boxes within a cell of the patch have a face on its nodes
    for (auto const& patch : level)
    {
        auto const near = core::grow(amr::phare_box_from<dimension>(patch->getBox()), 1);
        std::vector<core::Box<int, dimension>> nearby;
        for (auto const& box : covered)
            if (box * near)
                nearby.push_back(box);
        add(*patch, nearby);
    }
    return skipped;
}



/**
 * @brief BorderNode is an edge node of a refined level (regionEdgeNodes) with faces of E outside
 * the level: per direction d, bit 0 of outside[d] the face of E_d below it along d, bit 1 the
 * face above, `count` of them in all (1 on a side, 2 on a 2D outer corner, none on an inner one).
 */
template<std::size_t dimension>
struct BorderNode
{
    core::Point<int, dimension> node;
    std::array<std::uint8_t, dimension> outside{};
    std::uint8_t count = 0;
};

template<std::size_t dimension>
using PatchBorderNodes
    = std::unordered_map<std::string /*patch id*/, std::vector<BorderNode<dimension>>>;

/**
 * @brief levelBorderNodes returns the nodes of `cells` (AMR cell box) on the edge of `region`
 * (the level around the patch) that have faces of E outside it, a face being outside when none of
 * its 2^(dim-1) cells is in the region.
 */
template<typename Box_t>
auto levelBorderNodes(Box_t const& cells, std::vector<Box_t> const& region)
{
    constexpr auto dimension = Box_t::dimension;
    constexpr auto corners   = std::size_t{1} << dimension;

    // the face of E_d below (side 0) or above (side 1) the node: cell node_d - 1 + side along d,
    // node_e - 1 or node_e along the others
    auto const outside = [&](auto const& node, std::size_t const d, int const side) {
        for (std::size_t corner = 0; corner < corners; ++corner)
            if (((corner >> d) & 1u) == 0)
            {
                auto cell = node;
                cell[d] += side - 1;
                for (std::size_t e = 0; e < dimension; ++e)
                    if (e != d)
                        cell[e] -= static_cast<int>((corner >> e) & 1u);
                if (core::isIn(cell, region))
                    return false;
            }
        return true;
    };

    std::vector<BorderNode<dimension>> nodes;
    std::set<std::array<int, dimension>> seen; // regionEdgeNodes lists a node once per face
    for (auto const& box : regionEdgeNodes(cells, region))
    {
        BorderNode<dimension> border{box.lower};
        for (std::size_t d = 0; d < dimension; ++d)
            for (int const side : {0, 1})
                if (outside(border.node, d, side))
                {
                    border.outside[d] |= static_cast<std::uint8_t>(1u << side);
                    ++border.count;
                }
        if (border.count > 0 and seen.insert(border.node.toArray()).second)
            nodes.push_back(border);
    }
    return nodes;
}

/**
 * @brief levelBorderNodesOf returns, for every local patch of a refined level ln, its level border
 * nodes (levelBorderNodes) in its box and in its ghost box but its last ghost cell layer: those
 * whose faces it holds, as the particles of its band read them. None on the root.
 */
template<typename GridLayout>
auto levelBorderNodesOf(SAMRAI::hier::PatchHierarchy const& hierarchy, int const ln)
{
    constexpr auto dimension = GridLayout::dimension;
    constexpr auto ghosts    = static_cast<int>(GridLayout::options.field_ghost_width);

    PatchBorderNodes<dimension> nodes;
    if (ln == 0)
        return nodes;
    // the same-level boxes within the ghost box, periodic images included
    SAMRAI::hier::HierarchyNeighbors const neighbors{hierarchy, ln, ln, true, ghosts};
    for (auto const& patch : *hierarchy.getPatchLevel(ln))
    {
        auto const box = amr::phare_box_from<dimension>(patch->getBox());
        std::vector<core::Box<int, dimension>> region{box};
        for (auto const& neighbor : neighbors.getSameLevelNeighbors(patch->getBox(), ln))
            region.push_back(amr::phare_box_from<dimension>(neighbor));
        if (auto border = levelBorderNodes(core::grow(box, ghosts - 1), region); !border.empty())
            nodes[amr::to_string(patch->getGlobalId())] = std::move(border);
    }
    return nodes;
}

/**
 * @brief gaussOnLevelBorder sets the faces of E outside a refined level around each of its edge
 * nodes a patch holds (levelBorderNodesOf) so that Gauss's law holds on the node, div E = c2 rho,
 * with the model rho (depositFKCharge) at E's time: the residual of the node is taken off its n
 * outside faces, 1/n each. The fills give these faces the coarser E, while the charge of the edge
 * nodes follows the fine particles: their residual then changed at the rate of the coarse/fine
 * mismatch, and the cleaning mask kept it (BORDER-GAUSS, border_gauss.md). Set after every fill
 * from values every patch holding the node shares (the fills make the copies of E and rho equal
 * to their owner's), the faces are the same in every patch and hold no state: the residual of an
 * edge node is that of its charge noise only, at every step.
 */
template<typename GridLayout, typename Field, typename VecField>
void gaussOnLevelBorder(GridLayout const& layout, double const c2, Field const& rho, VecField& E,
                        std::vector<BorderNode<GridLayout::dimension>> const& nodes)
{
    constexpr auto dimension = GridLayout::dimension;
    auto const lower         = layout.AMRBox().lower;
    auto const dl            = layout.meshSize();

    for (auto const& border : nodes)
    {
        auto const& a = border.node;
        // the face of E_d below (side 0) or above (side 1) the node: dual along d
        auto const face = [&](std::size_t const d, int const side) {
            std::array<std::size_t, dimension> index;
            for (std::size_t e = 0; e < dimension; ++e)
                index[e] = layout.physicalStartIndex(E[d], static_cast<core::Direction>(e)) + a[e]
                           - lower[e] - (e == d ? 1 - side : 0);
            return index;
        };
        std::array<std::size_t, dimension> node;
        for (std::size_t d = 0; d < dimension; ++d)
            node[d]
                = layout.physicalStartIndex(rho, static_cast<core::Direction>(d)) + a[d] - lower[d];

        auto residual = -c2 * rho(node);
        for (std::size_t d = 0; d < dimension; ++d)
            residual += (E[d](face(d, 1)) - E[d](face(d, 0))) / dl[d];
        auto const share = residual / border.count;
        for (std::size_t d = 0; d < dimension; ++d)
            for (int const side : {0, 1})
                if ((border.outside[d] >> side) & 1u)
                    E[d](face(d, side)) += (side == 0 ? 1. : -1.) * share * dl[d];
    }
}

/**
 * @brief gaussOnLevelBorderOf is gaussOnLevelBorder on every patch of a level with border nodes
 */
template<typename FKModel, typename VecField>
void gaussOnLevelBorderOf(typename FKModel::level_t& level, FKModel& model, VecField& E,
                          PatchBorderNodes<FKModel::dimension> const& border)
{
    using GridLayout = FKModel::gridlayout_type;
    if (border.empty())
        return;
    auto& rho = model.state.rho;
    for (auto& patch : model.resourcesManager->enumerate(level, rho, E))
        if (auto const it = border.find(amr::to_string(patch->getGlobalId())); it != border.end())
            gaussOnLevelBorder(amr::layoutFromPatch<GridLayout>(*patch), model.c2(), rho, E,
                               it->second);
}



/**
 * @brief advanceFKElectric is steps 1 and 2 of SolverFK: J and rho from the particles at time
 * (depositFKCurrent), then E, at time - dt/2, cleaned (Marder, d from model.divergence_cleaning on
 * this level's mesh) and advanced to time + dt/2 (Ampere-Maxwell). E is the model E, or the copy of
 * it the solver advances one step ahead for a finer level (P4 plan D1). E ghosts are read, not
 * filled; labels are the caller's. The cleaning of a patch skips its nodes in `skipped`
 * (cleaningSkippedNodes; empty without boundary_cleaning_mask). On a refined level the faces of E
 * outside it first take Gauss's law on its edge nodes, `border` (gaussOnLevelBorder).
 */
template<typename FKModel, typename Messenger, typename VecField>
void advanceFKElectric(typename FKModel::level_t& level, FKModel& model, Messenger& messenger,
                       VecField& E, double const dt, double const time,
                       PatchBoxes<FKModel::dimension> const& skipped,
                       PatchBorderNodes<FKModel::dimension> const& border)
{
    using GridLayout = FKModel::gridlayout_type;

    auto& rm  = *model.resourcesManager;
    auto& B   = model.state.electromag.B;
    auto& J   = model.state.J;
    auto& rho = model.state.rho;

    depositFKCurrent(level, model, messenger, dt, /*back=*/1., /*front=*/1., time);
    gaussOnLevelBorderOf(level, model, E, border);

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
