

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
