#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/models/mhd_state.hpp"

#include "initializer/data_provider.hpp"

#include "tests/core/data/field/test_usable_field_fixtures_mhd.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <memory>
#include <vector>

using namespace PHARE::core;
using namespace PHARE::initializer;

// Vector-potential init of the MHD magnetic field: B = curl(A) with a full 3D vector potential
// A = (Ax, Ay, Az), built with the discrete `deriv` operator (the same Faraday uses).
//
// B = curl(A) is divergence-free in the discrete (Yee) sense BY CONSTRUCTION: div . curl = 0 for
// the staggered operators. Sampling the same analytic field at the face centres (component-wise
// init) is generically NOT discretely divergence-free.
//
// A *tilted* (non-separable) potential is used on purpose: for a grid-aligned separable mode the
// component-sampled init is accidentally discretely div-free, so it would not distinguish the two
// inits. With the coupled phase (kx x + ky y [+ kz z]) the B components are sampled at different
// staggered points and the component init has a genuine nonzero discrete divergence, while the
// discrete curl stays exactly div-free.
namespace
{
auto makeSpan(std::vector<double> vals)
{
    return std::make_shared<VectorSpan<double>>(std::move(vals));
}

template<std::size_t d>
VecFieldInitializer<d> makeAInit(InitFunction<d> ax, InitFunction<d> ay, InitFunction<d> az)
{
    PHAREDict dict;
    dict["x_component"] = ax;
    dict["y_component"] = ay;
    dict["z_component"] = az;
    return VecFieldInitializer<d>{dict};
}

InitFunction<2> const zero2 = [](auto const& x, auto const&) {
    return makeSpan(std::vector<double>(x.size(), 0.0));
};

// ---------------------------------------------------------------------------- 2D
constexpr std::size_t interp_order = 1;
constexpr std::uint32_t nx = 20, ny = 20;
constexpr double dx = 0.1, dy = 0.15;
// Distinct kx, ky (and dx != dy) on purpose: when kx==ky and dx==dy the component-sampled discrete
// divergence of a single plane wave cancels by symmetry and would look div-free too.
double const kx = 2.0 * M_PI / (nx * dx);
double const ky = 4.0 * M_PI / (ny * dy);

// out-of-plane vector potential A_z = psi = sin(kx x + ky y)
InitFunction<2> const psiInit = [](std::vector<double> const& x, std::vector<double> const& y) {
    std::vector<double> v(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        v[i] = std::sin(kx * x[i] + ky * y[i]);
    return makeSpan(std::move(v));
};

// analytic components of B = curl(psi z_hat): Bx = ky cos(.), By = -kx cos(.)
InitFunction<2> const bxInit = [](std::vector<double> const& x, std::vector<double> const& y) {
    std::vector<double> v(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        v[i] = ky * std::cos(kx * x[i] + ky * y[i]);
    return makeSpan(std::move(v));
};
InitFunction<2> const byInit = [](std::vector<double> const& x, std::vector<double> const& y) {
    std::vector<double> v(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        v[i] = -kx * std::cos(kx * x[i] + ky * y[i]);
    return makeSpan(std::move(v));
};

class MHDVectorPotentialInit : public ::testing::Test
{
protected:
    using GridLayout_t = GridLayout<GridLayoutImplYeeMHD<2, interp_order>>;

    GridLayout_t layout{{{dx, dy}}, {{nx, ny}}, Point{0.0, 0.0}};

    UsableVecFieldMHD<2> B{"B", layout, MHDQuantity::Vector::B};
    UsableVecFieldMHD<2> E{"E", layout, MHDQuantity::Vector::E};
    UsableFieldMHD<2> divB{"divB", layout, MHDQuantity::Scalar::divB};

    // max |div B| over the physical (cell-centred) box, with the staggered stencil
    // div B = dx Bx + dy By.
    double maxDivB(UsableVecFieldMHD<2>& Bvf)
    {
        auto const& Bx = Bvf(Component::X);
        auto const& By = Bvf(Component::Y);
        double m       = 0.0;
        layout.evalOnBox(divB.super(), [&](auto&... ijk) {
            double const d = layout.template deriv<Direction::X>(Bx, {ijk...})
                             + layout.template deriv<Direction::Y>(By, {ijk...});
            m = std::max(m, std::abs(d));
        });
        return m;
    }
};

} // namespace


// Legacy out-of-plane potential (A = (0, 0, psi)): reproduces the original 2D behavior exactly.
TEST_F(MHDVectorPotentialInit, curlInitIsDiscretelyDivergenceFree)
{
    auto aInit = makeAInit<2>(zero2, zero2, psiInit);
    initBFromPotential(aInit, B.super(), E.super(), layout);
    EXPECT_LT(maxDivB(B), 1e-10);
}

TEST_F(MHDVectorPotentialInit, componentSampledInitIsNotDivergenceFree)
{
    FieldUserFunctionInitializer::initialize(B(Component::X), layout, bxInit);
    FieldUserFunctionInitializer::initialize(B(Component::Y), layout, byInit);
    FieldUserFunctionInitializer::initialize(B(Component::Z), layout, zero2);
    EXPECT_GT(maxDivB(B), 1e-4);
}

TEST_F(MHDVectorPotentialInit, curlInitIsOrdersOfMagnitudeMoreDivFreeThanComponent)
{
    UsableVecFieldMHD<2> Bcomp{"Bcomp", layout, MHDQuantity::Vector::B};
    auto aInit = makeAInit<2>(zero2, zero2, psiInit);
    initBFromPotential(aInit, B.super(), E.super(), layout);
    FieldUserFunctionInitializer::initialize(Bcomp(Component::X), layout, bxInit);
    FieldUserFunctionInitializer::initialize(Bcomp(Component::Y), layout, byInit);
    FieldUserFunctionInitializer::initialize(Bcomp(Component::Z), layout, zero2);

    EXPECT_LT(maxDivB(B), 1e-6 * maxDivB(Bcomp));
}

// In-plane potential (A = (Ax, Ay, 0)) exercises the new Bz = dAy/dx - dAx/dy term: with only Az
// (the legacy path) Bz was forced to 0. Bz here must be a nonzero field.
TEST_F(MHDVectorPotentialInit, inPlanePotentialProducesNonzeroBz)
{
    auto aInit = makeAInit<2>(psiInit, psiInit, zero2);
    initBFromPotential(aInit, B.super(), E.super(), layout);

    auto const& Bz = B(Component::Z);
    double maxBz   = 0.0;
    layout.evalOnBox(divB.super(),
                     [&](auto&... ijk) { maxBz = std::max(maxBz, std::abs(Bz(ijk...))); });
    EXPECT_GT(maxBz, 1e-3); // ~ |kx - ky|
    // B = (0, 0, Bz) has no in-plane divergence contribution -> still discretely div-free.
    EXPECT_LT(maxDivB(B), 1e-10);
}


namespace
{
// ---------------------------------------------------------------------------- 3D
constexpr std::uint32_t nz = 20;
constexpr double dz            = 0.2;
double const kz                = 6.0 * M_PI / (nz * dz);

// Full 3D tilted vector potential, all three components nonzero, coupled phase.
InitFunction<3> const a3 = [](std::vector<double> const& x, std::vector<double> const& y,
                              std::vector<double> const& z) {
    std::vector<double> v(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        v[i] = std::sin(kx * x[i] + ky * y[i] + kz * z[i]);
    return makeSpan(std::move(v));
};

class MHDVectorPotentialInit3D : public ::testing::Test
{
protected:
    using GridLayout_t = GridLayout<GridLayoutImplYeeMHD<3, interp_order>>;

    GridLayout_t layout{{{dx, dy, dz}}, {{nx, ny, nz}}, Point{0.0, 0.0, 0.0}};

    UsableVecFieldMHD<3> B{"B", layout, MHDQuantity::Vector::B};
    UsableVecFieldMHD<3> E{"E", layout, MHDQuantity::Vector::E};
    UsableFieldMHD<3> divB{"divB", layout, MHDQuantity::Scalar::divB};

    // div B = dx Bx + dy By + dz Bz.
    double maxDivB(UsableVecFieldMHD<3>& Bvf)
    {
        auto const& Bx = Bvf(Component::X);
        auto const& By = Bvf(Component::Y);
        auto const& Bz = Bvf(Component::Z);
        double m       = 0.0;
        layout.evalOnBox(divB.super(), [&](auto&... ijk) {
            double const d = layout.template deriv<Direction::X>(Bx, {ijk...})
                             + layout.template deriv<Direction::Y>(By, {ijk...})
                             + layout.template deriv<Direction::Z>(Bz, {ijk...});
            m = std::max(m, std::abs(d));
        });
        return m;
    }
};

} // namespace


// Full 3D curl B = curl(A) with all three components active is discretely divergence-free.
TEST_F(MHDVectorPotentialInit3D, fullCurlInitIsDiscretelyDivergenceFree)
{
    auto aInit = makeAInit<3>(a3, a3, a3);
    initBFromPotential(aInit, B.super(), E.super(), layout);

    // sanity: the field is actually populated (nonzero), not trivially div-free by being all zero.
    auto const& Bx = B(Component::X);
    double maxB    = 0.0;
    layout.evalOnBox(divB.super(),
                     [&](auto&... ijk) { maxB = std::max(maxB, std::abs(Bx(ijk...))); });
    EXPECT_GT(maxB, 1e-3);

    EXPECT_LT(maxDivB(B), 1e-10);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
