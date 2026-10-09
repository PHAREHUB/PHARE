#include <cmath>
#include <cstdint>
#include <map>
#include <tuple>
#include <vector>
#include <sstream>
#include <algorithm>

#include "phare_mpi.hpp"

#include "core/utilities/types.hpp"
#include "core/data/particles/particle.hpp"
#include "amr/data/particles/refine/split.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


namespace
{
template<std::size_t dimension, std::size_t interpOrder, std::size_t refineParticlesNbr>
using Splitter
    = PHARE::amr::Splitter<PHARE::core::DimConst<dimension>, PHARE::core::InterpConst<interpOrder>,
                           PHARE::core::RefinedParticlesConst<refineParticlesNbr>>;

template<typename Splitter>
struct SplitterTest : public ::testing::Test
{
    SplitterTest() { Splitter splitter; }
};

// every hybrid permutation of res/sim/all.txt
using Splitters
    = testing::Types<Splitter<1, 1, 2>, Splitter<1, 1, 3>, Splitter<1, 2, 2>, Splitter<1, 2, 3>,
                     Splitter<1, 2, 4>, Splitter<1, 3, 2>, Splitter<1, 3, 3>, Splitter<1, 3, 4>,
                     Splitter<1, 3, 5>,                                                          //
                     Splitter<2, 1, 4>, Splitter<2, 1, 5>, Splitter<2, 1, 8>, Splitter<2, 1, 9>, //
                     Splitter<2, 2, 4>, Splitter<2, 2, 5>, Splitter<2, 2, 8>, Splitter<2, 2, 9>,
                     Splitter<2, 2, 16>, //
                     Splitter<2, 3, 4>, Splitter<2, 3, 5>, Splitter<2, 3, 8>, Splitter<2, 3, 9>,
                     Splitter<2, 3, 25>, //
                     Splitter<3, 1, 6>, Splitter<3, 1, 12>, Splitter<3, 1, 27>, Splitter<3, 2, 6>,
                     Splitter<3, 2, 12>, Splitter<3, 3, 6>, Splitter<3, 3, 12>>;
// Splitter<3, 2, 27> and Splitter<3, 3, 27> are disabled in res/sim/all.txt until their weights are
// known, see the TODO(@rochSmets) in split_3d.hpp

TYPED_TEST_SUITE(SplitterTest, Splitters);

TYPED_TEST(SplitterTest, constexpr_init)
{
    constexpr TypeParam param{};
}


// Smets et al. 2021, section 4 ("w0 + 4w1 = 1"): the children weights of one parent sum to 1
// (the dispatcher then multiplies by refinementRatio^dim, see splitter.hpp)
TYPED_TEST(SplitterTest, weights_sum_to_one)
{
    TypeParam splitter{};
    double sum = 0;
    PHARE::core::apply(splitter.patterns, [&](auto const& pattern) {
        sum += static_cast<double>(pattern.weight_) * pattern.deltas_.size();
    });
    EXPECT_NEAR(sum, 1., 1e-5);
}


// ParticlesRefineOperator::getSplitBox grows the destination box by splitBoxGrowth fine cells:
// any parent that puts a child in the box must lie within that distance of it
TYPED_TEST(SplitterTest, split_box_holds_every_contributing_parent)
{
    constexpr auto dim    = TypeParam::dimension;
    constexpr int growth  = PHARE::amr::splitBoxGrowth<TypeParam>();
    constexpr int maxCell = 4;

    TypeParam splitter{};
    std::array<PHARE::core::Particle<dim>, TypeParam::nbRefinedPart> children;
    int furthestContributor = 0;

    for (std::size_t iDim = 0; iDim < dim; ++iDim)
        for (int cell = -maxCell; cell <= maxCell; ++cell)
            for (double const delta : {0., 0.25, 0.5, 0.75, 0.999})
            {
                PHARE::core::Particle<dim> parent;
                parent.weight = 1;
                parent.delta.fill(0.5);
                parent.iCell[iDim] = cell;
                parent.delta[iDim] = delta;

                splitter(parent, children);
                for (auto const& child : children)
                    if (child.iCell[iDim] == 0) // child in the box, a single cell at the origin
                        furthestContributor = std::max(furthestContributor, std::abs(cell));
            }

    EXPECT_LE(furthestContributor, growth);
}


// (tau, w, delta) triplet of a subset of children, as defined in Smets et al. 2021, eq. (10):
// tau = a^2 + b^2 + c^2 for a child at (a delta, b delta, c delta), a, b, c in {-1, 0, 1}
struct Triplet
{
    int tau;
    double weight;
    double delta;
    std::size_t count;
};

template<typename Splitter_t>
std::vector<Triplet> triplets(Splitter_t const& splitter)
{
    std::vector<Triplet> out;
    PHARE::core::apply(splitter.patterns, [&](auto const& pattern) {
        auto const weight = static_cast<double>(pattern.weight_);
        for (auto const& d : pattern.deltas_)
        {
            int tau      = 0;
            double delta = 0;
            for (auto const x : d)
                if (std::abs(x) > 1e-6f)
                {
                    ++tau;
                    delta = std::max(delta, static_cast<double>(std::abs(x)));
                }
            auto it = std::find_if(out.begin(), out.end(), [&](auto const& t) {
                return t.tau == tau and std::abs(t.delta - delta) < 1e-5
                       and std::abs(t.weight - weight) < 1e-7;
            });
            if (it == out.end())
                out.push_back({tau, weight, delta, 1});
            else
                ++it->count;
        }
    });
    return out;
}

using Key = std::tuple<std::size_t, std::size_t, std::size_t>; // dim, interp, nbrRefinedPart

// Best-pattern (tau, w, delta) values of Smets et al. 2021, "A new method to dispatch split
// particles in Particle-In-Cell codes", CPC, arXiv:2104.10675.
// d=1: Tables 1 (w) and 2 (delta); d=2: Tables 3 and 4 (non-dagger rows); d=3: Tables 5 and 6.
// delta is on the refined grid. tau=0 children have delta 0.
std::map<Key, std::vector<Triplet>> const paperTriplets{
    // d = 1, Tables 1-2
    {{1, 1, 2}, {{1, 0.5, 0.551569, 2}}},
    {{1, 1, 3}, {{0, 0.5, 0, 1}, {1, 0.25, 1.0, 2}}},
    {{1, 2, 2}, {{1, 0.5, 0.663959, 2}}},
    {{1, 2, 3}, {{0, 0.468137, 0, 1}, {1, 0.265931, 1.112033, 2}}},
    {{1, 2, 4}, {{1, 0.375, 0.5, 2}, {1, 0.125, 1.5, 2}}},
    {{1, 3, 2}, {{1, 0.5, 0.752399, 2}}},
    {{1, 3, 3}, {{0, 0.473943, 0, 1}, {1, 0.263028, 1.275922, 2}}},
    {{1, 3, 4}, {{1, 0.364766, 0.542949, 2}, {1, 0.135234, 1.664886, 2}}},
    {{1, 3, 5}, {{0, 0.375, 0, 1}, {1, 0.25, 1.0, 2}, {1, 0.0625, 2.0, 2}}},
    // d = 2, Tables 3-4
    {{2, 1, 4}, {{2, 0.25, 0.571783, 4}}},
    {{2, 1, 5}, {{0, 0.239863, 0, 1}, {2, 0.190034, 0.721835, 4}}},
    {{2, 1, 8}, {{1, 0.179488, 0.700909, 4}, {2, 0.070512, 1.05786, 4}}},
    {{2, 1, 9}, {{0, 0.25, 0, 1}, {1, 0.125, 1.0, 4}, {2, 0.0625, 1.0, 4}}},
    {{2, 2, 4}, {{2, 0.25, 0.683734, 4}}},
    {{2, 2, 5}, {{0, 0.239166, 0, 1}, {1, 0.190209, 1.203227, 4}}},
    {{2, 2, 8}, {{1, 0.178624, 0.828428, 4}, {2, 0.071376, 1.236701, 4}}},
    {{2, 2, 9}, {{0, 0.213636, 0, 1}, {1, 0.126689, 1.105332, 4}, {2, 0.069902, 1.143884, 4}}},
    {{2, 3, 4}, {{2, 0.25, 0.776459, 4}}},
    {{2, 3, 5}, {{0, 0.242666, 0, 1}, {1, 0.189333, 1.376953, 4}}},
    {{2, 3, 8}, {{1, 0.179318, 0.942365, 4}, {2, 0.070682, 1.423324, 4}}},
    {{2, 3, 9}, {{0, 0.218605, 0, 1}, {1, 0.126871, 1.267689, 4}, {2, 0.068477, 1.315944, 4}}},
    // d = 3, Tables 5-6
    {{3, 1, 6}, {{1, 0.166666, 0.966431, 6}}},
    {{3, 2, 6}, {{1, 0.166666, 1.149658, 6}}},
    {{3, 3, 6}, {{1, 0.166666, 1.312622, 6}}},
    {{3, 1, 12}, {{2, 0.083333, 0.74823, 12}}},
    {{3, 2, 12}, {{2, 0.083333, 0.888184, 12}}},
    {{3, 3, 12}, {{2, 0.083333, 1.012756, 12}}},
    // d = 3, p = 1, N = 27 is the exact split ((p+2)^3 children, section 5): tensor product of
    // the d = 1 exact N = 3 split of Table 1 (w = 0.5, 0.25; delta = 1).
    // Tables 5-6 of arXiv:2104.10675 only list the tau = 0 and 1 rows for N = 27.
    {{3, 1, 27},
     {{0, 0.125, 0, 1}, {1, 0.0625, 1.0, 6}, {2, 0.03125, 1.0, 12}, {3, 0.015625, 1.0, 8}}},
};

TYPED_TEST(SplitterTest, matches_smets_2021_tables)
{
    Key const key{TypeParam::dimension, TypeParam::interp_order, TypeParam::nbRefinedPart};
    if (!paperTriplets.count(key))
        GTEST_SKIP() << "no reference values in Smets et al. 2021 for this configuration";

    auto const expected = paperTriplets.at(key);
    auto const actual   = triplets(TypeParam{});

    auto const describe = [](auto const& ts) {
        std::ostringstream os;
        for (auto const& t : ts)
            os << "{tau=" << t.tau << " w=" << t.weight << " delta=" << t.delta << " n=" << t.count
               << "} ";
        return os.str();
    };

    ASSERT_EQ(actual.size(), expected.size())
        << "expected " << describe(expected) << "\nactual   " << describe(actual);

    for (auto const& e : expected)
    {
        auto it = std::find_if(actual.begin(), actual.end(), [&](auto const& a) {
            return a.tau == e.tau and a.count == e.count and std::abs(a.weight - e.weight) < 1e-6
                   and std::abs(a.delta - e.delta) < 1e-5;
        });
        EXPECT_TRUE(it != actual.end()) << "missing {tau=" << e.tau << " w=" << e.weight
                                        << " delta=" << e.delta << " n=" << e.count << "}\n"
                                        << "actual   " << describe(actual);
    }
}

} // namespace
