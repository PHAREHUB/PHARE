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
#include <SAMRAI/hier/PatchHierarchy.h>

#include <algorithm>
#include <array>
#include <string>
#include <unordered_map>
#include <vector>


namespace PHARE::solver
{
// per direction d, faces of E_d (or J_d) as AMR indexes: the cell along d, the node along the
// others
template<std::size_t dimension>
using Faces = std::array<std::vector<core::Point<int, dimension>>, dimension>;

template<std::size_t dimension>
using PatchFaces = std::unordered_map<std::string /*patch id*/, Faces<dimension>>;

/**
 * @brief levelBorderFaces returns the faces of E_d a patch holds just outside `region` (the patch
 * and its same-level neighbours, nonLevelGhostBox), d the direction across the side: the first
 * ghost face of each side, where no cell along it (2^(dim-1) around the face) is in the region,
 * corners included. Gauss's law on the edge nodes of a refined level reads them, and no patch of
 * the level advances them. Refilled from the coarser E, they changed the residual of those nodes
 * at the rate of the coarse/fine mismatch, which the cleaning mask kept (BORDER-GAUSS,
 * border_gauss.md): the solver advances them with the level's Ampere (ampereOnFaces) and keeps
 * them through the ghost fills, so that the residual of an edge node changes as its charge does.
 */
template<typename Box_t>
auto levelBorderFaces(Box_t const& cells, std::vector<Box_t> const& region)
{
    constexpr auto dimension = Box_t::dimension;
    constexpr auto corners   = std::size_t{1} << dimension;

    Faces<dimension> faces;
    for (std::size_t d = 0; d < dimension; ++d)
        for (auto const cell : {cells.lower[d] - 1, cells.upper[d] + 1})
        {
            auto side = cells; // the nodes along the other directions
            for (std::size_t e = 0; e < dimension; ++e)
                side.upper[e] += e != d;
            side.lower[d] = side.upper[d] = cell;
            core::for_each_point(side, [&](auto const& face) {
                bool outside = true;
                for (std::size_t corner = 0; corner < corners; ++corner)
                    if (((corner >> d) & 1u) == 0)
                    {
                        auto around = face;
                        for (std::size_t e = 0; e < dimension; ++e)
                            around[e] -= static_cast<int>((corner >> e) & 1u);
                        outside = outside and !core::isIn(around, region);
                    }
                if (outside)
                    faces[d].push_back(face);
            });
        }
    return faces;
}

// the local index in `field` (E_d or J_d) of the face of AMR index `face`
template<typename GridLayout, typename Field, typename Point_t>
auto localFace(GridLayout const& layout, Field const& field, Point_t const& face)
{
    std::array<std::size_t, GridLayout::dimension> index;
    for (std::size_t d = 0; d < GridLayout::dimension; ++d)
        index[d] = layout.physicalStartIndex(field, static_cast<core::Direction>(d)) + face[d]
                   - layout.AMRBox().lower[d];
    return index;
}

// one binomial pass (1/4, 1/2, 1/4 in every direction) of `field` at the local `index`, as
// core::BinomialFilter computes it on the physical nodes
template<typename Field, std::size_t dimension>
double binomialAt(Field const& field, std::array<std::size_t, dimension> const& index)
{
    constexpr std::array<double, 3> w{0.25, 0.5, 0.25};
    std::size_t stencil = 1;
    for (std::size_t d = 0; d < dimension; ++d)
        stencil *= 3;
    double sum = 0.;
    for (std::size_t k = 0; k < stencil; ++k)
    {
        auto at       = index;
        double weight = 1.;
        for (std::size_t d = 0, rest = k; d < dimension; ++d, rest /= 3)
        {
            at[d] = at[d] + rest % 3 - 1;
            weight *= w[rest % 3];
        }
        sum += weight * field(at);
    }
    return sum;
}

/**
 * @brief ampereOnFaces is the Ampere-Maxwell update of AmpereMaxwell on the given faces of E (the
 * level border faces of the patch), E_d += dt c2 ((curl B)_d - J_d): their B ghosts and J (band
 * deposit) exist there.
 */
template<typename GridLayout, typename VecField>
void ampereOnFaces(GridLayout const& layout, double const c2, VecField const& B, VecField const& J,
                   VecField& E, double const dt, Faces<GridLayout::dimension> const& faces)
{
    auto const curlB = core::curl(layout, B);
    auto const coef  = dt * c2;
    for (std::size_t d = 0; d < GridLayout::dimension; ++d)
        for (auto const& face : faces[d])
        {
            auto const index = localFace(layout, E[d], face);
            auto const mesh  = core::MeshIndex<GridLayout::dimension>{index};
            auto const curl  = d == 0 ? curlB.x(mesh) : d == 1 ? curlB.y(mesh) : curlB.z(mesh);
            E[d](index) += coef * (curl - J[d](index));
        }
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
                      double const dt, double const back, double const front, double const time,
                      PatchFaces<FKModel::dimension> const& borderFaces = {})
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
    // on its own nodes, so a node shared by patches gets it once from each (P4 plan D5). Each
    // patch adds these terms in its own order: the max makes the copies equal again
    if (level.getLevelNumber() > 0)
    {
        deposit([](auto const& pop) -> auto const& { return pop.levelGhostParticles(); });
        messenger.equalizeCurrentBorders(J, level, time);
        if (cleaning)
            messenger.equalizeChargeDensityBorders(rho, level, rhoTime);
    }

    // fill then filter, every pass
    auto const filter = [&](auto& field, auto&& fillGhosts) {
        for (std::size_t pass = 0; pass < model.current_filter_passes; ++pass)
        {
            fillGhosts();
            for (auto& patch : rm.enumerate(level, field))
                core::BinomialFilter{amr::layoutFromPatch<GridLayout>(*patch)}(field);
        }
    };

    // J as `filter`, and its level border faces (levelBorderFaces) filtered as physical faces,
    // from the band deposit around them: the residual of an edge node reads the filtered charge
    // and the J its outside face is advanced with (ampereOnFaces). Exact at N = 1; at N > 1 the
    // second ghost layer is not filtered again, as the band nodes
    for (std::size_t pass = 0; pass < model.current_filter_passes; ++pass)
    {
        messenger.fillCurrentGhosts(J, level, time);
        for (auto& patch : rm.enumerate(level, J))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            auto const it     = borderFaces.find(amr::to_string(patch->getGlobalId()));
            std::array<std::vector<double>, FKModel::dimension> faceJ;
            if (it != borderFaces.end())
                for (std::size_t d = 0; d < FKModel::dimension; ++d)
                    for (auto const& face : it->second[d])
                        faceJ[d].push_back(binomialAt(J[d], localFace(layout, J[d], face)));
            core::BinomialFilter{layout}(J);
            if (it != borderFaces.end())
                for (std::size_t d = 0; d < FKModel::dimension; ++d)
                    for (std::size_t f = 0; f < faceJ[d].size(); ++f)
                        J[d](localFace(layout, J[d], it->second[d][f])) = faceJ[d][f];
        }
    }
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
using PatchBoxes
    = std::unordered_map<std::string /*patch id*/, std::vector<core::Box<int, dimension>>>;


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
 * @brief levelBorderFacesOf returns levelBorderFaces for every local patch of a refined level ln
 * that has any; none on the root.
 */
template<typename GridLayout>
auto levelBorderFacesOf(SAMRAI::hier::PatchHierarchy const& hierarchy, int const ln)
{
    constexpr auto dimension = GridLayout::dimension;
    PatchFaces<dimension> faces;
    if (ln == 0)
        return faces;
    auto const neighbors = amr::makeSameLevelNeighbors(hierarchy, ln);
    for (auto const& patch : *hierarchy.getPatchLevel(ln))
    {
        auto patchFaces
            = levelBorderFaces(amr::phare_box_from<dimension>(patch->getBox()),
                               amr::makeNonLevelGhostBoxFor<GridLayout>(*patch, neighbors));
        if (std::any_of(patchFaces.begin(), patchFaces.end(),
                        [](auto const& v) { return !v.empty(); }))
            faces[amr::to_string(patch->getGlobalId())] = std::move(patchFaces);
    }
    return faces;
}



/**
 * @brief advanceFKElectric is steps 1 and 2 of SolverFK: J and rho from the particles at time
 * (depositFKCurrent), then E, at time - dt/2, cleaned (Marder, d from model.divergence_cleaning on
 * this level's mesh) and advanced to time + dt/2 (Ampere-Maxwell). E is the model E, or the copy of
 * it the solver advances one step ahead for a finer level (P4 plan D1). E ghosts are read, not
 * filled; labels are the caller's. The cleaning of a patch skips its nodes in `skipped`
 * (cleaningSkippedNodes; empty without boundary_cleaning_mask). On a refined level the Ampere
 * update also advances the level border faces of E (`borderFaces`, levelBorderFaces), and J is
 * filtered there too.
 */
template<typename FKModel, typename Messenger, typename VecField>
void advanceFKElectric(typename FKModel::level_t& level, FKModel& model, Messenger& messenger,
                       VecField& E, double const dt, double const time,
                       PatchBoxes<FKModel::dimension> const& skipped,
                       PatchFaces<FKModel::dimension> const& borderFaces)
{
    using GridLayout = FKModel::gridlayout_type;

    auto& rm  = *model.resourcesManager;
    auto& B   = model.state.electromag.B;
    auto& J   = model.state.J;
    auto& rho = model.state.rho;

    depositFKCurrent(level, model, messenger, dt, /*back=*/1., /*front=*/1., time, borderFaces);

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
        if (auto const it = borderFaces.find(amr::to_string(patch->getGlobalId()));
            it != borderFaces.end())
            ampereOnFaces(layout, model.c2(), B, J, E, dt, it->second);
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
