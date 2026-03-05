
#include "test_linear_combinaisons_yee.hpp"

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
    if (dim == Dim && order == Order)
    {
        using Layout = GridLayout<GridLayoutImplYee<Dim, Order>>;
        auto actual  = call.template operator()<Layout>();
        compareCombi<Dim>(combi, actual);
    }
}

template<typename Callable>
void runTestFile(std::string const& filename, Callable&& call)
{
    auto expectedCombinaisons = readFile(filename);

    for (auto const& combi : expectedCombinaisons)
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

TEST(MomentToEx, combinaisonOk)
{
    runTestFile("linear_coefs_yee_momentToEx.txt",
                []<typename Layout>() { return Layout::momentsToEx(); });
}

TEST(MomentToEy, combinaisonOk)
{
    runTestFile("linear_coefs_yee_momentToEy.txt",
                []<typename Layout>() { return Layout::momentsToEy(); });
}

TEST(MomentToEz, combinaisonOk)
{
    runTestFile("linear_coefs_yee_momentToEz.txt",
                []<typename Layout>() { return Layout::momentsToEz(); });
}

TEST(ExToMoment, combinaisonOk)
{
    runTestFile("linear_coefs_yee_ExToMoment.txt",
                []<typename Layout>() { return Layout::ExToMoments(); });
}

TEST(EyToMoment, combinaisonOk)
{
    runTestFile("linear_coefs_yee_EyToMoment.txt",
                []<typename Layout>() { return Layout::EyToMoments(); });
}

TEST(EzToMoment, combinaisonOk)
{
    runTestFile("linear_coefs_yee_EzToMoment.txt",
                []<typename Layout>() { return Layout::EzToMoments(); });
}

TEST(JxToMoment, combinaisonOk)
{
    runTestFile("linear_coefs_yee_JxToMoment.txt",
                []<typename Layout>() { return Layout::JxToMoments(); });
}

TEST(JyToMoment, combinaisonOk)
{
    runTestFile("linear_coefs_yee_JyToMoment.txt",
                []<typename Layout>() { return Layout::JyToMoments(); });
}

TEST(JzToMoment, combinaisonOk)
{
    runTestFile("linear_coefs_yee_JzToMoment.txt",
                []<typename Layout>() { return Layout::JzToMoments(); });
}

TEST(ByToEx, combinaisonOk)
{
    runTestFile("linear_coefs_yee_ByToEx.txt", []<typename Layout>() { return Layout::ByToEx(); });
}

TEST(BzToEx, combinaisonOk)
{
    runTestFile("linear_coefs_yee_BzToEx.txt", []<typename Layout>() { return Layout::BzToEx(); });
}

TEST(BxToEy, combinaisonOk)
{
    runTestFile("linear_coefs_yee_BxToEy.txt", []<typename Layout>() { return Layout::BxToEy(); });
}

TEST(BzToEy, combinaisonOk)
{
    runTestFile("linear_coefs_yee_BzToEy.txt", []<typename Layout>() { return Layout::BzToEy(); });
}

TEST(BxToEz, combinaisonOk)
{
    runTestFile("linear_coefs_yee_BxToEz.txt", []<typename Layout>() { return Layout::BxToEz(); });
}

TEST(ByToEz, combinaisonOk)
{
    runTestFile("linear_coefs_yee_ByToEz.txt", []<typename Layout>() { return Layout::ByToEz(); });
}

TEST(JxToEx, combinaisonOk)
{
    runTestFile("linear_coefs_yee_JxToEx.txt", []<typename Layout>() { return Layout::JxToEx(); });
}

TEST(JyToEy, combinaisonOk)
{
    runTestFile("linear_coefs_yee_JyToEy.txt", []<typename Layout>() { return Layout::JyToEy(); });
}

TEST(JzToEz, combinaisonOk)
{
    runTestFile("linear_coefs_yee_JzToEz.txt", []<typename Layout>() { return Layout::JzToEz(); });
}

TEST(BxToEx, combinaisonOk)
{
    runTestFile("linear_coefs_yee_BxToEx.txt", []<typename Layout>() { return Layout::BxToEx(); });
}

TEST(ByToEy, combinaisonOk)
{
    runTestFile("linear_coefs_yee_ByToEy.txt", []<typename Layout>() { return Layout::ByToEy(); });
}

TEST(BzToEz, combinaisonOk)
{
    runTestFile("linear_coefs_yee_BzToEz.txt", []<typename Layout>() { return Layout::BzToEz(); });
}
