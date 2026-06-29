#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/models/mhd_state.hpp"

#include "tests/core/data/field/test_usable_field_fixtures_mhd.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <memory>
#include <vector>

using namespace PHARE::core;
using namespace PHARE::initializer;

// Vector-potential init of the MHD magnetic field (2D).
//
// B = curl(A_z z_hat) built with the discrete `deriv` operator (the same Faraday uses) is
// divergence-free in the discrete (Yee) sense BY CONSTRUCTION: div . curl = 0 for the staggered
// operators. Sampling the same analytic field at the face centres (component-wise init) is
// generically NOT discretely divergence-free.
//
// A *tilted* (non-separable) potential is used on purpose: for a grid-aligned separable mode
// like cos(kx x) cos(ky y) the component-sampled init is accidentally discretely div-free, so it
// would not distinguish the two inits. With the coupled phase (kx x + ky y) the two B components
// are sampled at different staggered points and the component init has a genuine nonzero discrete
// divergence, while the discrete curl stays exactly div-free.
//
// psi(x,y) = sin(kx x + ky y)  =>  B = curl(psi z_hat) = (dpsi/dy, -dpsi/dx, 0):
//   Bx =  dpsi/dy =  ky cos(kx x + ky y)
//   By = -dpsi/dx = -kx cos(kx x + ky y)
namespace
{
constexpr std::size_t dim          = 2;
constexpr std::size_t interp_order = 1;
constexpr std::uint32_t nx = 20, ny = 20;
constexpr double dx = 0.1, dy = 0.15;
// Distinct kx, ky (and dx != dy) on purpose: when kx==ky and dx==dy the component-sampled
// discrete divergence of a single plane wave cancels by symmetry and would look div-free too.
double const kx = 2.0 * M_PI / (nx * dx);
double const ky = 4.0 * M_PI / (ny * dy);

auto makeSpan(std::vector<double> vals)
{
    return std::make_shared<VectorSpan<double>>(std::move(vals));
}

// out-of-plane vector potential A_z = psi
InitFunction<dim> const psiInit = [](std::vector<double> const& x, std::vector<double> const& y) {
    std::vector<double> v(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        v[i] = std::sin(kx * x[i] + ky * y[i]);
    return makeSpan(std::move(v));
};

// analytic components of B = curl(psi z_hat)
InitFunction<dim> const bxInit = [](std::vector<double> const& x, std::vector<double> const& y) {
    std::vector<double> v(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        v[i] = ky * std::cos(kx * x[i] + ky * y[i]);
    return makeSpan(std::move(v));
};
InitFunction<dim> const byInit = [](std::vector<double> const& x, std::vector<double> const& y) {
    std::vector<double> v(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
        v[i] = -kx * std::cos(kx * x[i] + ky * y[i]);
    return makeSpan(std::move(v));
};
InitFunction<dim> const bzInit = [](std::vector<double> const& x, std::vector<double> const&) {
    return makeSpan(std::vector<double>(x.size(), 0.0));
};

class MHDVectorPotentialInit : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYeeMHD<dim, interp_order>;
    using GridLayout_t   = GridLayout<GridLayoutImpl>;

    GridLayout_t layout{{{dx, dy}}, {{nx, ny}}, Point{0.0, 0.0}};

    UsableVecFieldMHD<dim> B{"B", layout, MHDQuantity::Vector::B};
    UsableVecFieldMHD<dim> E{"E", layout, MHDQuantity::Vector::E};
    UsableFieldMHD<dim> divB{"divB", layout, MHDQuantity::Scalar::divB};

    // max |div B| over the physical (cell-centred) box, with the staggered stencil
    // div B = dx Bx + dy By.
    template<typename VecField>
    double maxDivB(VecField& Bvf)
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


TEST_F(MHDVectorPotentialInit, curlInitIsDiscretelyDivergenceFree)
{
    initBFromPotentialAz(psiInit, B.super(), E(Component::Z), layout);
    EXPECT_LT(maxDivB(B), 1e-10);
}

TEST_F(MHDVectorPotentialInit, componentSampledInitIsNotDivergenceFree)
{
    FieldUserFunctionInitializer::initialize(B(Component::X), layout, bxInit);
    FieldUserFunctionInitializer::initialize(B(Component::Y), layout, byInit);
    FieldUserFunctionInitializer::initialize(B(Component::Z), layout, bzInit);
    EXPECT_GT(maxDivB(B), 1e-4);
}

TEST_F(MHDVectorPotentialInit, curlInitIsOrdersOfMagnitudeMoreDivFreeThanComponent)
{
    UsableVecFieldMHD<dim> Bcomp{"Bcomp", layout, MHDQuantity::Vector::B};
    initBFromPotentialAz(psiInit, B.super(), E(Component::Z), layout);
    FieldUserFunctionInitializer::initialize(Bcomp(Component::X), layout, bxInit);
    FieldUserFunctionInitializer::initialize(Bcomp(Component::Y), layout, byInit);
    FieldUserFunctionInitializer::initialize(Bcomp(Component::Z), layout, bzInit);

    EXPECT_LT(maxDivB(B), 1e-6 * maxDivB(Bcomp));
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
