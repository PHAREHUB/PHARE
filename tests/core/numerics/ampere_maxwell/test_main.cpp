#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/numerics/ampere/ampere_maxwell.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <numbers>


using namespace PHARE::core;


template<std::size_t dim_, std::size_t interp_order_>
struct Opts
{
    static constexpr std::size_t dim          = dim_;
    static constexpr std::size_t interp_order = interp_order_;
};


template<typename Opts_t>
class AmpereMaxwellTest : public ::testing::Test
{
protected:
    static constexpr std::size_t dim          = Opts_t::dim;
    static constexpr std::size_t interp_order = Opts_t::interp_order;
    using GridLayout_t = PHARE_Types<PHARE::SimOpts{dim, interp_order}>::Hybrid::GridLayout_t;

    static constexpr std::uint32_t nCells = 32;
    static constexpr double dx            = 0.1;
    static constexpr double L             = nCells * dx;
    static constexpr double c2            = 25.;
    static constexpr double dt            = 0.01;
    static constexpr double k             = 2 * std::numbers::pi / L;
    static constexpr double sentinel      = 12345.;

    GridLayout_t layout{ConstArray<double, dim>(dx), ConstArray<std::uint32_t, dim>(nCells),
                        Point<double, dim>{}};
    UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
    UsableVecField<dim> J{"J", layout, HybridQuantity::Vector::J};
    UsableVecField<dim> E{"E", layout, HybridQuantity::Vector::E};


    auto position(auto const& f, auto const... ijk) const
    {
        return layout.fieldNodeCoordinates(
            f, layout.localToAMR(Point{static_cast<std::uint32_t>(ijk)...}.as_signed()));
    }

    bool isPhysical(auto const& f, auto const... ijk) const
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

    void fillAll(auto& vecfield, double value)
    {
        for (auto& f : vecfield)
            layout.evalOnGhostBox(f, [&](auto const&... ijk) { f(ijk...) = value; });
    }


    // B_c = sum over the simulated directions d of sin(k x_d + phase(c, d)), so that every
    // derivative entering the curl is non zero
    static double phase(std::size_t c, std::size_t d) { return 0.3 * (c + 1) + 0.7 * d; }

    static double bValue(std::size_t c, auto const& r)
    {
        double value = 0.;
        for (std::size_t d = 0; d < dim; ++d)
            value += std::sin(k * r[d] + phase(c, d));
        return value;
    }

    static double bDerivative(std::size_t c, std::size_t d, auto const& r)
    {
        return d < dim ? k * std::cos(k * r[d] + phase(c, d)) : 0.;
    }

    // (curl B)_c = d_i B_j - d_j B_i, (c, i, j) a circular permutation of (x, y, z)
    static double curlB(std::size_t c, auto const& r)
    {
        auto const i = (c + 1) % 3;
        auto const j = (c + 2) % 3;
        return bDerivative(j, i, r) - bDerivative(i, j, r);
    }
};

using AllOpts = ::testing::Types<Opts<1, 1>, Opts<1, 2>, Opts<1, 3>, Opts<2, 1>, Opts<2, 2>,
                                 Opts<2, 3>, Opts<3, 1>, Opts<3, 2>, Opts<3, 3>>;
TYPED_TEST_SUITE(AmpereMaxwellTest, AllOpts);


TYPED_TEST(AmpereMaxwellTest, currentOnlyTermIsExact)
{
    this->fillAll(this->B, 0.);
    for (auto c : {Component::X, Component::Y, Component::Z})
    {
        auto& Jc = this->J(c);
        this->layout.evalOnGhostBox(
            Jc, [&](auto const&... ijk) { Jc(ijk...) = 1. + static_cast<double>(c); });
    }
    this->fillAll(this->E, 0.);

    AmpereMaxwell{this->layout, this->c2}(this->B, this->J, this->E, this->dt);

    for (auto c : {Component::X, Component::Y, Component::Z})
    {
        auto& Ec = this->E(c);
        this->layout.evalOnBox(Ec, [&](auto const&... ijk) {
            EXPECT_DOUBLE_EQ(Ec(ijk...), -this->dt * this->c2 * (1. + static_cast<double>(c)));
        });
    }
}


TYPED_TEST(AmpereMaxwellTest, curlTermIsSecondOrderAndGhostsUntouched)
{
    this->fillAll(this->J, 0.);
    this->fillAll(this->E, this->sentinel);
    for (auto& f : this->E)
        this->layout.evalOnBox(f, [&](auto const&... ijk) { f(ijk...) = 0.; });

    for (std::size_t c = 0; c < 3; ++c)
    {
        auto& Bc = this->B[c];
        this->layout.evalOnGhostBox(Bc, [&](auto const&... ijk) {
            Bc(ijk...) = this->bValue(c, this->position(Bc, ijk...));
        });
    }

    AmpereMaxwell{this->layout, this->c2}(this->B, this->J, this->E, this->dt);

    // each Yee derivative of sin(k x) is centered: error k^3 dx^2 / 24, two per component
    auto const amp = this->dt * this->c2 * this->k;
    auto const kdx = this->k * this->dx;
    auto const tol = amp * kdx * kdx / 6;

    for (std::size_t c = 0; c < 3; ++c)
    {
        auto& Ec = this->E[c];
        this->layout.evalOnGhostBox(Ec, [&](auto const&... ijk) {
            if (this->isPhysical(Ec, ijk...))
            {
                auto const expected
                    = this->dt * this->c2 * this->curlB(c, this->position(Ec, ijk...));
                EXPECT_NEAR(Ec(ijk...), expected, tol) << "component " << c;
            }
            else
                EXPECT_EQ(Ec(ijk...), this->sentinel) << "component " << c;
        });
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
