#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/numerics/ampere/ampere_maxwell.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <numbers>


using namespace PHARE::core;


template<typename InterpOrder>
class AmpereMaxwell1D : public ::testing::Test
{
protected:
    static constexpr std::size_t dim          = 1;
    static constexpr std::size_t interp_order = InterpOrder::value;
    using GridLayout_t = PHARE_Types<PHARE::SimOpts{dim, interp_order}>::Hybrid::GridLayout_t;

    static constexpr std::size_t nCells = 100;
    static constexpr double dx          = 0.1;
    static constexpr double L           = nCells * dx;
    static constexpr double c2          = 25.;
    static constexpr double dt          = 0.01;
    static constexpr double k           = 2 * std::numbers::pi / L;
    static constexpr double sentinel    = 12345.;

    GridLayout_t layout{{{dx}}, {{nCells}}, Point{0.}};
    UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
    UsableVecField<dim> J{"J", layout, HybridQuantity::Vector::J};
    UsableVecField<dim> E{"E", layout, HybridQuantity::Vector::E};

    // physical coordinate of local node i of field f
    double x(auto const& f, std::uint32_t i) const
    {
        auto const start   = layout.physicalStartIndex(f, Direction::X);
        auto const isDual  = layout.centering(f.physicalQuantity())[0] == QtyCentering::dual;
        auto const shifted = static_cast<double>(i) - static_cast<double>(start);
        return (shifted + (isDual ? 0.5 : 0.)) * dx;
    }

    void fillAll(auto& vecfield, double value)
    {
        for (auto& f : vecfield)
            layout.evalOnGhostBox(f, [&](auto const&... ijk) { f(ijk...) = value; });
    }

    bool isPhysical(auto const& f, std::uint32_t i) const
    {
        return i >= layout.physicalStartIndex(f, Direction::X)
               and i <= layout.physicalEndIndex(f, Direction::X);
    }
};

using InterpOrders = ::testing::Types<std::integral_constant<std::size_t, 1>,
                                      std::integral_constant<std::size_t, 2>,
                                      std::integral_constant<std::size_t, 3>>;
TYPED_TEST_SUITE(AmpereMaxwell1D, InterpOrders);


TYPED_TEST(AmpereMaxwell1D, currentOnlyTermIsExact)
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


TYPED_TEST(AmpereMaxwell1D, curlTermIsSecondOrderAndGhostsUntouched)
{
    auto const k = this->k;
    this->fillAll(this->J, 0.);
    this->fillAll(this->B, 0.);
    this->fillAll(this->E, this->sentinel);
    for (auto& f : this->E)
        this->layout.evalOnBox(f, [&](auto const&... ijk) { f(ijk...) = 0.; });

    auto& By = this->B(Component::Y);
    auto& Bz = this->B(Component::Z);
    this->layout.evalOnGhostBox(By, [&](auto const& i) { By(i) = std::sin(k * this->x(By, i)); });
    this->layout.evalOnGhostBox(Bz, [&](auto const& i) { Bz(i) = std::cos(k * this->x(Bz, i)); });

    AmpereMaxwell{this->layout, this->c2}(this->B, this->J, this->E, this->dt);

    // (curl B)_y = -dBz/dx = k sin(kx), (curl B)_z = dBy/dx = k cos(kx)
    auto const amp = this->dt * this->c2 * k;
    auto const tol = amp * (k * this->dx) * (k * this->dx); // second order finite difference

    auto& Ex = this->E(Component::X);
    auto& Ey = this->E(Component::Y);
    auto& Ez = this->E(Component::Z);

    this->layout.evalOnGhostBox(Ey, [&](auto const& i) {
        if (this->isPhysical(Ey, i))
            EXPECT_NEAR(Ey(i), amp * std::sin(k * this->x(Ey, i)), tol);
        else
            EXPECT_EQ(Ey(i), this->sentinel);
    });
    this->layout.evalOnGhostBox(Ez, [&](auto const& i) {
        if (this->isPhysical(Ez, i))
            EXPECT_NEAR(Ez(i), amp * std::cos(k * this->x(Ez, i)), tol);
        else
            EXPECT_EQ(Ez(i), this->sentinel);
    });
    this->layout.evalOnGhostBox(Ex, [&](auto const& i) {
        if (this->isPhysical(Ex, i))
            EXPECT_EQ(Ex(i), 0.);
        else
            EXPECT_EQ(Ex(i), this->sentinel);
    });
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
