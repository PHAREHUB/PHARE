// BorderSums adds the contributions to a node in the order of the contributing patch boxes, then
// patch ids, from zero: the result is a function of the set of contributions, the same on every
// patch holding the node, whatever the order the transactions delivered them in.
//
// Four 2D patches meet at a corner, so that nodes around it are held by four patches. Each patch
// stores the overlaps of its ghost box with the others, in some registration order, then folds.
// Every ghost box node must equal its contributions added in that order from zero.

#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/models/options/hybrid_options.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/utilities/array.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/utilities/point/point.hpp"
#include "phare_core.hpp"

#include "phare_mpi.hpp" // IWYU pragma: keep

#include "amr/data/field/field_data.hpp"
#include "amr/data/field/field_overlap.hpp"
#include "amr/messengers/field_border_sum.hpp"

#include "simulator/simulator_def.hpp"

#include <SAMRAI/hier/BlockId.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/hier/GlobalId.h>
#include <SAMRAI/hier/Index.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/LocalId.h>
#include <SAMRAI/hier/Transformation.h>
#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <optional>
#include <random>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

using namespace PHARE::core;
using namespace PHARE::amr;

namespace
{
std::size_t constexpr dim    = 2;
std::size_t constexpr interp = 1;

using GridLayout_t = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;
using Grid_t       = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
using FieldData_t  = FieldData<GridLayout_t, Grid_t>;
using Box_t        = PHARE::core::Box<int, dim>;

auto constexpr qty     = HybridQuantity::Scalar::rho; // primal, like the moments summed
int constexpr dataId   = 0;
int constexpr cellSize = 4;

struct Patch
{
    Box_t box;
    SAMRAI::hier::Box samraiBox;
    std::shared_ptr<FieldData_t> data;

    Box_t ghostNodes() const { return data->gridLayout.AMRGhostBoxFor(qty); }
};

SAMRAI::hier::Box toSamrai(Box_t const& box, int const localId = 0)
{
    SAMRAI::tbox::Dimension const d{dim};
    SAMRAI::hier::Index lower{d}, upper{d};
    for (std::size_t i = 0; i < dim; ++i)
    {
        lower[i] = box.lower[i];
        upper[i] = box.upper[i];
    }
    return {lower, upper, SAMRAI::hier::BlockId{0}, SAMRAI::hier::LocalId{localId},
            SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank()};
}

Box_t cellBox(int const ix, int const iy)
{
    return {{ix * cellSize, iy * cellSize},
            {ix * cellSize + cellSize - 1, iy * cellSize + cellSize - 1}};
}

// four patches around the corner at cell (cellSize, cellSize)
std::vector<Box_t> cornerBoxes()
{
    return {cellBox(0, 0), cellBox(1, 0), cellBox(0, 1), cellBox(1, 1)};
}

// two distinct patches with one box, and a neighbour
std::vector<Box_t> identicalBoxes()
{
    return {cellBox(0, 0), cellBox(0, 0), cellBox(1, 0)};
}

// one patch per box, with patch id its index, values spread over many orders of magnitude so that
// the order of the additions shows in the result
std::vector<Patch> makePatches(std::vector<Box_t> const& boxes)
{
    std::mt19937_64 gen{2026};
    std::uniform_real_distribution<double> mantissa{1., 2.};
    std::uniform_int_distribution<int> exponent{-30, 30};

    SAMRAI::tbox::Dimension const d{dim};
    std::vector<Patch> patches;
    int localId = 0;
    for (auto const& box : boxes)
    {
        GridLayout_t const layout{ConstArray<double, dim>(0.1),
                                  ConstArray<std::uint32_t, dim>(cellSize),
                                  Point<double, dim>{0.1 * box.lower[0], 0.1 * box.lower[1]}, box};
        auto const samraiBox = toSamrai(box, localId++);
        auto data            = std::make_shared<FieldData_t>(
            samraiBox, SAMRAI::hier::IntVector{d, GridLayout_t::options.field_ghost_width}, "rho",
            layout, qty);
        for_each_point(layout.AMRGhostBoxFor(qty), [&](auto const& node) {
            data->field(layout.AMRToLocal(node))
                = std::ldexp(mantissa(gen), exponent(gen)) * (gen() % 2 ? 1. : -1.);
        });
        patches.push_back({box, samraiBox, data});
    }
    return patches;
}

double valueAt(Patch const& patch, auto const& node)
{
    return patch.data->field(patch.data->gridLayout.AMRToLocal(node));
}

struct Transaction
{
    std::size_t dst, src;
    std::shared_ptr<FieldOverlap> overlap;
};

std::vector<Transaction> makeTransactions(std::vector<Patch> const& patches)
{
    SAMRAI::tbox::Dimension const d{dim};
    std::vector<Transaction> transactions;
    for (std::size_t dst = 0; dst < patches.size(); ++dst)
        for (std::size_t src = 0; src < patches.size(); ++src)
            if (dst != src)
                if (auto const overlap = patches[dst].ghostNodes() * patches[src].ghostNodes())
                    transactions.push_back(
                        {dst, src,
                         std::make_shared<FieldOverlap>(
                             SAMRAI::hier::BoxContainer{toSamrai(*overlap)},
                             SAMRAI::hier::Transformation{SAMRAI::hier::IntVector::getZero(d)})});
    return transactions;
}

// the patches of `boxes` after the border sum, transactions registered in the given order
std::vector<Patch> borderSum(std::vector<Box_t> const& boxes,
                             std::vector<Transaction> const& transactions)
{
    auto patches = makePatches(boxes);
    BorderSums<FieldData_t> sums;

    std::vector<std::size_t> handles;
    for (auto const& t : transactions)
        handles.push_back(
            sums.add(patches[t.dst].samraiBox, patches[t.src].samraiBox, t.overlap, dataId));

    auto const patchOf = [&](auto const& id) -> Patch& {
        return *std::ranges::find_if(
            patches, [&](auto const& p) { return p.samraiBox.getGlobalId() == id; });
    };
    sums.finalize([&](auto const& id) -> auto const& { return patchOf(id).samraiBox; });

    // as the schedule would: copies of the sources before any patch is summed
    for (std::size_t i = 0; i < transactions.size(); ++i)
        patches[transactions[i].src].data->packTo(sums.contribution(handles[i]),
                                                  *transactions[i].overlap);

    sums.fold([&](auto const& id, int) -> FieldData_t& { return *patchOf(id).data; });

    return patches;
}

// every ghost box node of every patch equals its contributions added from zero in the order of the
// contributing patch boxes, then patch ids; returns the number of nodes held by all patches
std::size_t expectSumsInSourceOrder(std::vector<Box_t> const& boxes)
{
    auto const original = makePatches(boxes);
    auto const summed   = borderSum(boxes, makeTransactions(original));

    auto const sourceLess = [](Patch const* a, Patch const* b) {
        if (!(a->box == b->box))
            return lexicographicLess(a->box, b->box);
        return a->samraiBox.getGlobalId() < b->samraiBox.getGlobalId();
    };

    std::size_t nodesHeldByAll = 0;
    for (std::size_t ip = 0; ip < original.size(); ++ip)
    {
        auto const& patch = original[ip];
        for_each_point(patch.ghostNodes(), [&](auto const& node) {
            std::vector<Patch const*> holders;
            for (auto const& other : original)
                if (isIn(node, other.ghostNodes()))
                    holders.push_back(&other);
            std::ranges::sort(holders, sourceLess);
            nodesHeldByAll += holders.size() == original.size();

            double expected = 0;
            for (auto const* holder : holders)
                expected += valueAt(*holder, node);

            auto const got = valueAt(summed[ip], node);
            EXPECT_EQ(std::memcmp(&got, &expected, sizeof(double)), 0)
                << "patch " << ip << " node " << node << ": " << got << " vs " << expected;
        });
    }
    return nodesHeldByAll;
}

// the same bits on every patch whatever the order the transactions are registered in
void expectIndependentOfTransactionOrder(std::vector<Box_t> const& boxes)
{
    auto transactions    = makeTransactions(makePatches(boxes));
    auto const reference = borderSum(boxes, transactions);

    std::mt19937 gen{7};
    for (int shuffle = 0; shuffle < 5; ++shuffle)
    {
        std::ranges::shuffle(transactions, gen);
        auto const summed = borderSum(boxes, transactions);
        for (std::size_t ip = 0; ip < reference.size(); ++ip)
        {
            auto const& field = reference[ip].data->field;
            EXPECT_EQ(std::memcmp(summed[ip].data->field.data(), field.data(),
                                  field.size() * sizeof(double)),
                      0)
                << "patch " << ip << ", shuffle " << shuffle;
        }
    }
}

} // namespace



TEST(BorderSums, addsTheContributionsToANodeInPatchBoxOrderFromZero)
{
    EXPECT_GT(expectSumsInSourceOrder(cornerBoxes()), 0u);
}


TEST(BorderSums, resultDoesNotDependOnTheTransactionOrder)
{
    expectIndependentOfTransactionOrder(cornerBoxes());
}


TEST(BorderSums, ordersIdenticalBoxesOfDistinctPatchesByPatchId)
{
    EXPECT_GT(expectSumsInSourceOrder(identicalBoxes()), 0u);
    expectIndependentOfTransactionOrder(identicalBoxes());

    // the two patches with one box hold the same nodes: they get the same bits
    auto const boxes  = identicalBoxes();
    auto const summed = borderSum(boxes, makeTransactions(makePatches(boxes)));
    auto const& field = summed[0].data->field;
    EXPECT_EQ(
        std::memcmp(summed[1].data->field.data(), field.data(), field.size() * sizeof(double)), 0);
}


TEST(BorderSums, throwsOnTwoContributionsFromOnePatch)
{
    auto patches            = makePatches(cornerBoxes());
    auto const transactions = makeTransactions(patches);
    auto const& t           = transactions.front();

    BorderSums<FieldData_t> sums;
    for (int i = 0; i < 2; ++i)
        sums.add(patches[t.dst].samraiBox, patches[t.src].samraiBox, t.overlap, dataId);

    EXPECT_THROW(
        sums.finalize([&](auto const&) -> auto const& { return patches[t.dst].samraiBox; }),
        std::runtime_error);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    int testResult = RUN_ALL_TESTS();

    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
