

#include "phare_core.hpp"
#include "phare_simulator_options.hpp"

#include "test_linear_combinations_yee.hpp"

#include "gtest/gtest.h"

#include <algorithm>

using namespace PHARE::core;

template<auto dimension>
struct Ordered
{
    std::array<int, dimension> idx;
    double coef;

    auto operator<=>(Ordered const&) const = default;
};

template<auto dimension>
auto makeExpected(auto const& combi)
{
    std::vector<Ordered<dimension>> result;

    for (std::size_t i = 0; i < combi.ix.size(); ++i)
    {
        Ordered<dimension> k{};

        if constexpr (dimension >= 1)
            k.idx[0] = combi.ix[i];
        if constexpr (dimension >= 2)
            k.idx[1] = combi.iy[i];
        if constexpr (dimension >= 3)
            k.idx[2] = combi.iz[i];

        k.coef = combi.coef;

        result.push_back(k);
    }

    return result;
}

template<auto dimension, typename ActualContainer>
auto makeActual(ActualContainer const& actual)
{
    std::vector<Ordered<dimension>> result;

    for (auto const& p : actual)
    {
        Ordered<dimension> k{};

        for (std::size_t d = 0; d < dimension; ++d)
            k.idx[d] = p.indexes[d];

        k.coef = p.coef;

        result.push_back(k);
    }

    return result;
}

template<auto dimension, typename ActualContainer>
void compareCombi(auto const& combi, ActualContainer const& actual)
{
    auto expected = makeExpected<dimension>(combi);
    auto obtained = makeActual<dimension>(actual);

    ASSERT_EQ(expected.size(), obtained.size());

    std::ranges::sort(expected);
    std::ranges::sort(obtained);

    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        for (std::size_t d = 0; d < dimension; ++d)
            EXPECT_EQ(expected[i].idx[d], obtained[i].idx[d]);

        EXPECT_DOUBLE_EQ(expected[i].coef, obtained[i].coef);
    }
}

template<int Dim, int Order, typename Callable>
void tryDispatch(int dim, int order, Callable&& call, auto const& combi)
{
    using GridLayoutT = PHARE_Types<PHARE::SimOpts{Dim, Order}>::Hybrid::GridLayout_t;
    if (dim == Dim && order == Order)
    {
        auto actual = call.template operator()<GridLayoutT>();
        compareCombi<Dim>(combi, actual);
    }
}

template<typename Callable>
void runTestFile(std::string const& filename, Callable&& call)
{
    auto expectedCombinations = readFile(filename);

    for (auto const& combi : expectedCombinations)
    {
        int dim   = combi.dimension;
        int order = combi.interpOrder;

        tryDispatch<1, 1>(dim, order, call, combi);
        tryDispatch<1, 2>(dim, order, call, combi);
        tryDispatch<1, 3>(dim, order, call, combi);

        tryDispatch<2, 1>(dim, order, call, combi);
        tryDispatch<2, 2>(dim, order, call, combi);
        tryDispatch<2, 3>(dim, order, call, combi);

        tryDispatch<3, 1>(dim, order, call, combi);
        tryDispatch<3, 2>(dim, order, call, combi);
        tryDispatch<3, 3>(dim, order, call, combi);
    }
}

TEST(MomentToEx, combinationOk)
{
    runTestFile("linear_coefs_yee_momentToEx.txt",
                []<typename Layout>() { return Layout::momentsToEx(); });
}

TEST(MomentToEy, combinationOk)
{
    runTestFile("linear_coefs_yee_momentToEy.txt",
                []<typename Layout>() { return Layout::momentsToEy(); });
}

TEST(MomentToEz, combinationOk)
{
    runTestFile("linear_coefs_yee_momentToEz.txt",
                []<typename Layout>() { return Layout::momentsToEz(); });
}

TEST(ExToMoment, combinationOk)
{
    runTestFile("linear_coefs_yee_ExToMoment.txt",
                []<typename Layout>() { return Layout::ExToMoments(); });
}

TEST(EyToMoment, combinationOk)
{
    runTestFile("linear_coefs_yee_EyToMoment.txt",
                []<typename Layout>() { return Layout::EyToMoments(); });
}

TEST(EzToMoment, combinationOk)
{
    runTestFile("linear_coefs_yee_EzToMoment.txt",
                []<typename Layout>() { return Layout::EzToMoments(); });
}

TEST(JxToMoment, combinationOk)
{
    runTestFile("linear_coefs_yee_JxToMoment.txt",
                []<typename Layout>() { return Layout::JxToMoments(); });
}

TEST(JyToMoment, combinationOk)
{
    runTestFile("linear_coefs_yee_JyToMoment.txt",
                []<typename Layout>() { return Layout::JyToMoments(); });
}

TEST(JzToMoment, combinationOk)
{
    runTestFile("linear_coefs_yee_JzToMoment.txt",
                []<typename Layout>() { return Layout::JzToMoments(); });
}

TEST(ByToEx, combinationOk)
{
    runTestFile("linear_coefs_yee_ByToEx.txt", []<typename Layout>() { return Layout::ByToEx(); });
}

TEST(BzToEx, combinationOk)
{
    runTestFile("linear_coefs_yee_BzToEx.txt", []<typename Layout>() { return Layout::BzToEx(); });
}

TEST(BxToEy, combinationOk)
{
    runTestFile("linear_coefs_yee_BxToEy.txt", []<typename Layout>() { return Layout::BxToEy(); });
}

TEST(BzToEy, combinationOk)
{
    runTestFile("linear_coefs_yee_BzToEy.txt", []<typename Layout>() { return Layout::BzToEy(); });
}

TEST(BxToEz, combinationOk)
{
    runTestFile("linear_coefs_yee_BxToEz.txt", []<typename Layout>() { return Layout::BxToEz(); });
}

TEST(ByToEz, combinationOk)
{
    runTestFile("linear_coefs_yee_ByToEz.txt", []<typename Layout>() { return Layout::ByToEz(); });
}

TEST(JxToEx, combinationOk)
{
    runTestFile("linear_coefs_yee_JxToEx.txt", []<typename Layout>() { return Layout::JxToEx(); });
}

TEST(JyToEy, combinationOk)
{
    runTestFile("linear_coefs_yee_JyToEy.txt", []<typename Layout>() { return Layout::JyToEy(); });
}

TEST(JzToEz, combinationOk)
{
    runTestFile("linear_coefs_yee_JzToEz.txt", []<typename Layout>() { return Layout::JzToEz(); });
}

TEST(BxToEx, combinationOk)
{
    runTestFile("linear_coefs_yee_BxToEx.txt", []<typename Layout>() { return Layout::BxToEx(); });
}

TEST(ByToEy, combinationOk)
{
    runTestFile("linear_coefs_yee_ByToEy.txt", []<typename Layout>() { return Layout::ByToEy(); });
}

TEST(BzToEz, combinationOk)
{
    runTestFile("linear_coefs_yee_BzToEz.txt", []<typename Layout>() { return Layout::BzToEz(); });
}


// ---------------------------------------------------------------------------
// Primitive B — dual-direction prolongation (coarse->fine, ratio 2).
// Self-contained algebraic checks (no Python-generated reference files).
// ---------------------------------------------------------------------------

using PHARE::core::dirX;
using ImplYee1 = PHARE_Types<PHARE::SimOpts{1, 1}>::Hybrid::GridLayout_t::implT;

namespace
{
    // sum of a stencil's coefficients
    auto coefSum(auto const& row)
    {
        double s = 0.;
        for (auto const& wp : row)
            s += wp.coef;
        return s;
    }

    // coefficient at a given offset along dirX (0 if absent)
    auto coefAt(auto const& row, int offset)
    {
        for (auto const& wp : row)
            if (wp.indexes[dirX] == offset)
                return wp.coef;
        return 0.;
    }

    // value the stencil produces from coarse cell-averages u(I+offset) = f(I+offset)
    auto applyRow(auto const& row, int I, auto&& f)
    {
        double v = 0.;
        for (auto const& wp : row)
            v += wp.coef * f(I + wp.indexes[dirX]);
        return v;
    }
} // namespace

// weights of each ladder rung match prolongation_operators.md §2.2 exactly
TEST(DualProlongation, order2WeightsMatchLadder)
{
    auto right = ImplYee1::directionalProlongation<dirX, +1, 2>(); // σ=+1
    EXPECT_DOUBLE_EQ(coefAt(right, -1), -1. / 8.);
    EXPECT_DOUBLE_EQ(coefAt(right, 0), 1.);
    EXPECT_DOUBLE_EQ(coefAt(right, +1), 1. / 8.);
}

// every rung is consistent: weights sum to 1 (reproduces a constant field)
TEST(DualProlongation, weightsSumToOne)
{
    EXPECT_DOUBLE_EQ(coefSum(ImplYee1::directionalProlongation<dirX, +1, 0>()), 1.);
    EXPECT_DOUBLE_EQ(coefSum(ImplYee1::directionalProlongation<dirX, -1, 0>()), 1.);
    EXPECT_DOUBLE_EQ(coefSum(ImplYee1::directionalProlongation<dirX, +1, 2>()), 1.);
    EXPECT_DOUBLE_EQ(coefSum(ImplYee1::directionalProlongation<dirX, -1, 2>()), 1.);
}

// conservation: the two children average back to ū_I at every order
TEST(DualProlongation, childrenMeanBackToCoarse)
{
    auto check = [](auto const& right, auto const& left) {
        // offset 0 carries the full coarse value in each child -> mean 1
        EXPECT_DOUBLE_EQ(0.5 * (coefAt(right, 0) + coefAt(left, 0)), 1.);
        // antisymmetric correction -> every off-center offset cancels in the mean
        for (int o : {-2, -1, 1, 2})
            EXPECT_DOUBLE_EQ(0.5 * (coefAt(right, o) + coefAt(left, o)), 0.);
    };
    check(ImplYee1::directionalProlongation<dirX, +1, 2>(),
          ImplYee1::directionalProlongation<dirX, -1, 2>());
}

// exactness: on linear cell-averages u_i = a + b*i, each child reproduces a + b*(I ± 1/4)
TEST(DualProlongation, exactOnLinearData)
{
    constexpr double a = 0.7, b = -1.3;
    auto f             = [](int i) { return a + b * i; };
    int const I        = 5;

    for (auto const& [row, sigma] :
         {std::pair{ImplYee1::directionalProlongation<dirX, +1, 2>(), +1},
          std::pair{ImplYee1::directionalProlongation<dirX, -1, 2>(), -1}})
        EXPECT_DOUBLE_EQ(applyRow(row, I, f), a + b * (I + sigma * 0.25));
}
