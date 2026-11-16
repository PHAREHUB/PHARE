#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/numerics/filter/binomial_filter.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <numbers>
#include <vector>


using namespace PHARE::core;


template<std::size_t dim_, std::size_t interp_order_>
struct Case
{
    static constexpr std::size_t dim          = dim_;
    static constexpr std::size_t interp_order = interp_order_;
};


// even numbers of cells, for the Nyquist mode, and mesh sizes that differ in every direction
template<typename Case_t>
class BinomialFilterTest : public ::testing::Test
{
protected:
    static constexpr std::size_t dim          = Case_t::dim;
    static constexpr std::size_t interp_order = Case_t::interp_order;
    using GridLayout_t = PHARE_Types<PHARE::SimOpts{dim, interp_order}>::Hybrid::GridLayout_t;

    static constexpr std::array<double, 3> allMeshSizes{0.1, 0.13, 0.07};
    static constexpr std::array<std::uint32_t, 3> allCells{16, 12, 10};

    static auto sub(auto const& all)
    {
        std::array<typename std::decay_t<decltype(all)>::value_type, dim> a;
        std::copy_n(all.begin(), dim, a.begin());
        return a;
    }

    GridLayout_t layout{sub(allMeshSizes), sub(allCells), Point<double, dim>{}};
    // J and B: every Yee centering, dual and primal in each direction
    UsableVecField<dim> J{"J", layout, HybridQuantity::Vector::J};
    UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};

    BinomialFilter<GridLayout_t> filter{layout};

    double length(std::size_t d) const { return allCells[d] * allMeshSizes[d]; }

    auto position(auto const& f, auto const&... ijk) const
    {
        return layout.fieldNodeCoordinates(
            f, layout.localToAMR(Point{static_cast<std::uint32_t>(ijk)...}.as_signed()));
    }

    bool isPhysical(auto const& f, auto const&... ijk) const
    {
        std::array<std::uint32_t, dim> const index{static_cast<std::uint32_t>(ijk)...};
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const direction = static_cast<Direction>(d);
            if (index[d] < layout.physicalStartIndex(f, direction)
                or index[d] > layout.physicalEndIndex(f, direction))
                return false;
        }
        return true;
    }

    // fn(position) on every node, ghosts included
    void fill(auto& f, auto&& fn) const
    {
        layout.evalOnGhostBox(f, [&](auto const&... ijk) { f(ijk...) = fn(position(f, ijk...)); });
    }

    auto copy(auto const& f) const { return std::vector<double>(f.data(), f.data() + f.size()); }

    // prod_d cos(k_d x_d + phase_d), k_d = 2 pi modes_d / L_d
    auto cosine(std::array<int, dim> const& modes) const
    {
        return [this, modes](auto const& r) {
            double value = 1.;
            for (std::size_t d = 0; d < dim; ++d)
                value *= std::cos(2 * std::numbers::pi * modes[d] / length(d) * r[d] + 0.3 + d);
            return value;
        };
    }

    // prod_d cos^2(k_d dl_d / 2), the transfer function of one pass
    double transfer(std::array<int, dim> const& modes) const
    {
        double t = 1.;
        for (std::size_t d = 0; d < dim; ++d)
            t *= std::pow(std::cos(std::numbers::pi * modes[d] / allCells[d]), 2);
        return t;
    }

    // one pass on a cosine: physical nodes multiplied by the transfer function, ghosts untouched
    void checkCosine(auto& f, std::array<int, dim> const& modes)
    {
        fill(f, cosine(modes));
        auto const before = copy(f);
        filter(f);
        checkFiltered(f, before, modes);
    }

    void checkFiltered(auto const& f, std::vector<double> const& before,
                       std::array<int, dim> const& modes) const
    {
        auto const t = transfer(modes);
        layout.evalOnGhostBox(f, [&](auto const&... ijk) {
            auto const i = &f(ijk...) - f.data();
            if (isPhysical(f, ijk...))
                EXPECT_NEAR(f(ijk...), t * before[i], 1e-14)
                    << f.name() << " modes " << modes[0] << " transfer " << t;
            else
                EXPECT_EQ(f(ijk...), before[i]) << f.name() << " ghost modified";
        });
    }
};

using Cases = ::testing::Types<Case<1, 1>, Case<1, 2>, Case<1, 3>, Case<2, 1>, Case<2, 2>,
                               Case<2, 3>, Case<3, 1>, Case<3, 2>, Case<3, 3>>;
TYPED_TEST_SUITE(BinomialFilterTest, Cases);



// along each direction alone: low modes and the Nyquist mode, which one pass removes
TYPED_TEST(BinomialFilterTest, transferAlongEachDirection)
{
    constexpr auto dim = TestFixture::dim;
    for (auto* vecfield : {&this->J, &this->B})
        for (auto& f : *vecfield)
            for (std::size_t d = 0; d < dim; ++d)
                for (int const mode : {1, 3, static_cast<int>(this->allCells[d] / 2)})
                {
                    std::array<int, dim> modes{};
                    modes[d] = mode;
                    this->checkCosine(f, modes);
                }
}


// along all directions at once: the stencil is the tensor product, edges and corners included
TYPED_TEST(BinomialFilterTest, transferAlongAllDirections)
{
    constexpr auto dim = TestFixture::dim;
    std::array<int, dim> modes;
    for (std::size_t d = 0; d < dim; ++d)
        modes[d] = static_cast<int>(d) + 1;

    for (auto* vecfield : {&this->J, &this->B})
        for (auto& f : *vecfield)
            this->checkCosine(f, modes);
}


TYPED_TEST(BinomialFilterTest, constantIsUnchanged)
{
    for (auto& f : this->J)
    {
        this->fill(f, [](auto const&) { return 1.7; });
        this->filter(f);
        this->layout.evalOnGhostBox(
            f, [&](auto const&... ijk) { EXPECT_DOUBLE_EQ(f(ijk...), 1.7) << f.name(); });
    }
}


// the vector field overload filters each component
TYPED_TEST(BinomialFilterTest, vecFieldFiltersEachComponent)
{
    constexpr auto dim = TestFixture::dim;
    std::array<std::array<int, dim>, 3> modes;
    std::array<std::vector<double>, 3> before;
    for (std::size_t c = 0; c < 3; ++c)
    {
        for (std::size_t d = 0; d < dim; ++d)
            modes[c][d] = static_cast<int>((c + d) % 3) + 1;
        this->fill(this->J[c], this->cosine(modes[c]));
        before[c] = this->copy(this->J[c]);
    }

    this->filter(this->J);

    for (std::size_t c = 0; c < 3; ++c)
        this->checkFiltered(this->J[c], before[c], modes[c]);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
