#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>

using namespace PHARE::core;

auto constexpr primal = QtyCentering::primal;
auto constexpr dual   = QtyCentering::dual;


template<std::size_t dim>
constexpr std::array<QtyCentering, dim> centeringFromBits(unsigned bits)
{
    std::array<QtyCentering, dim> c{};
    for (std::size_t d = 0; d < dim; ++d)
        c[d] = ((bits >> d) & 1u) ? dual : primal;
    return c;
}

template<std::size_t dim, typename Fn>
void forEachPair(Fn&& fn)
{
    constexpr unsigned nCenterings = 1u << dim;
    for_N<nCenterings * nCenterings>([&](auto ic) {
        constexpr unsigned i = ic();
        constexpr auto src   = centeringFromBits<dim>(i % nCenterings);
        constexpr auto dst   = centeringFromBits<dim>(i / nCenterings);
        fn.template operator()<src, dst>();
    });
}

template<typename Fn>
void forEachPairAllDims(Fn&& fn)
{
    forEachPair<1>(fn);
    forEachPair<2>(fn);
    forEachPair<3>(fn);
}

template<std::size_t dim>
std::array<int, dim> toArray(Point<int, dim> const& p)
{
    std::array<int, dim> a{};
    for (std::size_t d = 0; d < dim; ++d)
        a[d] = p[d];
    return a;
}


TEST(RecenteringStencil, LinearOffsets)
{
    EXPECT_EQ(linearRecenteringOffsets(primal), (std::array<int, 2>{0, 1}));
    EXPECT_EQ(linearRecenteringOffsets(dual), (std::array<int, 2>{-1, 0}));
}


TEST(RecenteringStencil, SameCenteringIsIdentity)
{
    forEachPairAllDims([]<auto src, auto dst>() {
        if constexpr (src == dst)
        {
            constexpr auto stencil = recenteringStencil<src, src>();
            EXPECT_EQ(stencil.size(), 1u);
            for (std::size_t d = 0; d < src.size(); ++d)
                EXPECT_EQ(stencil[0].indexes[d], 0);
            EXPECT_EQ(stencil[0].coef, 1.0);
        }
    });
}


TEST(RecenteringStencil, AllPairsShapeAndOrder)
{
    forEachPairAllDims([]<auto src, auto dst>() {
        constexpr std::size_t dim = src.size();
        constexpr auto stencil    = recenteringStencil<src, dst>();

        std::size_t p = 0;
        for (std::size_t d = 0; d < dim; ++d)
            p += (src[d] != dst[d]);

        EXPECT_EQ(stencil.size(), std::size_t{1} << p);

        double sum = 0.;
        for (auto const& wp : stencil)
        {
            sum += wp.coef;
            for (std::size_t d = 0; d < dim; ++d)
            {
                if (src[d] == dst[d])
                    EXPECT_EQ(wp.indexes[d], 0);
                else
                {
                    auto const offsets = linearRecenteringOffsets(src[d]);
                    EXPECT_TRUE(wp.indexes[d] == offsets[0] || wp.indexes[d] == offsets[1]);
                }
            }
        }
        EXPECT_EQ(sum, 1.0);

        for (std::size_t k = 1; k < stencil.size(); ++k)
        {
            auto const prev = toArray(stencil[k - 1].indexes);
            auto const cur  = toArray(stencil[k].indexes);
            EXPECT_TRUE(
                std::lexicographical_compare(prev.begin(), prev.end(), cur.begin(), cur.end()));
        }
    });
}


TEST(RecenteringStencil, ReproducesLinearFunctions)
{
    auto constexpr slopes = std::array{1.0, 0.5, 0.25};

    forEachPairAllDims([&]<auto src, auto dst>() {
        constexpr std::size_t dim = src.size();
        constexpr auto stencil    = recenteringStencil<src, dst>();
        int constexpr center      = 5;

        auto const position = [&](auto const& centering, std::array<int, dim> const& idx) {
            double f = 1.0;
            for (std::size_t d = 0; d < dim; ++d)
                f += slopes[d] * (idx[d] + (centering[d] == dual ? 0.5 : 0.0));
            return f;
        };

        double interpolated = 0.;
        for (auto const& wp : stencil)
        {
            std::array<int, dim> idx{};
            for (std::size_t d = 0; d < dim; ++d)
                idx[d] = center + wp.indexes[d];
            interpolated += wp.coef * position(src, idx);
        }

        std::array<int, dim> destIdx{};
        for (std::size_t d = 0; d < dim; ++d)
            destIdx[d] = center;
        EXPECT_DOUBLE_EQ(interpolated, position(dst, destIdx));
    });
}


TEST(RecenteringStencil, Order4ReproducesCubics)
{
    auto const cubic = [](double x) { return x * x * x - 2. * x * x + x; };

    auto const check = [&]<auto src, auto dst>() {
        constexpr auto stencil = recenteringStencil<src, dst, 4>();
        EXPECT_EQ(stencil.size(), 4u);

        int constexpr center = 5;
        auto const pos       = [&](QtyCentering c, int i) { return i + (c == dual ? 0.5 : 0.0); };

        double interpolated = 0.;
        for (auto const& wp : stencil)
            interpolated += wp.coef * cubic(pos(src[0], center + wp.indexes[0]));

        EXPECT_DOUBLE_EQ(interpolated, cubic(pos(dst[0], center)));
    };

    check.template operator()<std::array{primal}, std::array{dual}>();
    check.template operator()<std::array{dual}, std::array{primal}>();
}


TEST(RecenteringStencil, TensorProductOrder)
{
    constexpr auto a = std::array{WeightPoint<2>{Point<int, 2>{0, 0}, 0.5},
                                  WeightPoint<2>{Point<int, 2>{1, 0}, 0.5}};
    constexpr auto b = std::array{WeightPoint<2>{Point<int, 2>{0, -1}, 0.5},
                                  WeightPoint<2>{Point<int, 2>{0, 0}, 0.5}};

    constexpr auto result = tensorProduct<0, 1>(a, b);
    ASSERT_EQ(result.size(), 4u);

    constexpr std::array<std::array<int, 2>, 4> expected{{{0, -1}, {0, 0}, {1, -1}, {1, 0}}};
    for (std::size_t k = 0; k < 4; ++k)
    {
        EXPECT_EQ(result[k].indexes[0], expected[k][0]) << k;
        EXPECT_EQ(result[k].indexes[1], expected[k][1]) << k;
        EXPECT_EQ(result[k].coef, 0.25) << k;
    }
}


using YeeLayout2D = PHARE_Types<PHARE::SimOpts{2, 1}>::Hybrid::GridLayout_t;
using Grid2D      = Grid<NdArrayVector<2>, HybridQuantity::Scalar>;

TEST(RecenteringStencil, ProjectAppliesStencil)
{
    auto const layout = TestGridLayout<YeeLayout2D>::make(10);
    Grid2D bx{"Bx", layout, HybridQuantity::Scalar::Bx, 0.};
    auto& f = *bx;

    for (std::uint32_t i = 0; i < f.shape()[0]; ++i)
        for (std::uint32_t j = 0; j < f.shape()[1]; ++j)
            f(i, j) = 1.0 * i * i + 3.0 * j;

    constexpr auto bxCentering = YeeLayout2D::centering(HybridQuantity::Scalar::Bx);
    constexpr auto exCentering = YeeLayout2D::centering(HybridQuantity::Scalar::Ex);

    Point<std::size_t, 2> const index{5, 5};

    constexpr auto stencil = recenteringStencil<bxCentering, exCentering>();
    double expected        = 0.;
    for (auto const& wp : stencil)
        expected += wp.coef * f(index + wp.indexes);

    EXPECT_EQ((project<bxCentering, exCentering>(f, index)), expected);
}


struct FormerTable
{
    std::string model, name, src, dst;
    std::size_t dim, n;
    std::vector<double> values;
};

std::vector<FormerTable> readFormerTables()
{
    std::vector<FormerTable> tables;
    std::ifstream file(PHARE_RECENTERING_REFERENCE);
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream in(line);
        FormerTable t;
        in >> t.model >> t.name >> t.dim >> t.src >> t.dst >> t.n;
        double v;
        while (in >> v)
            t.values.push_back(v);
        tables.push_back(std::move(t));
    }
    return tables;
}

template<std::size_t dim>
std::string centeringString(std::array<QtyCentering, dim> const& c)
{
    std::string s;
    for (auto v : c)
        s += (v == primal ? 'p' : 'd');
    return s;
}


TEST(RecenteringStencil, MatchesFormerTablesInOrder)
{
    auto const tables = readFormerTables();
    ASSERT_EQ(tables.size(), 174u);

    for (auto const& table : tables)
    {
        std::size_t matched = 0;
        forEachPairAllDims([&]<auto src, auto dst>() {
            constexpr std::size_t dim = src.size();
            if (dim != table.dim || centeringString(src) != table.src
                || centeringString(dst) != table.dst)
                return;
            ++matched;

            constexpr auto stencil  = recenteringStencil<src, dst>();
            std::string const label = table.model + " " + table.name + " " + std::to_string(dim);
            ASSERT_EQ(stencil.size(), table.n) << label;
            ASSERT_EQ(table.values.size(), table.n * (dim + 1)) << label;
            for (std::size_t k = 0; k < stencil.size(); ++k)
            {
                for (std::size_t d = 0; d < dim; ++d)
                    EXPECT_EQ(stencil[k].indexes[d],
                              static_cast<int>(table.values[k * (dim + 1) + d]))
                        << label << " point " << k;
                EXPECT_EQ(stencil[k].coef, table.values[k * (dim + 1) + dim])
                    << label << " point " << k;
            }
        });
        EXPECT_EQ(matched, 1u) << table.model << " " << table.name;
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
