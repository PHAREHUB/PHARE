#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/numerics/differential_operators.hpp"

#include "gtest/gtest.h"

#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <stdexcept>

namespace PHARE::core
{

constexpr SimOpts mhdOpts(std::size_t dim, MHDOpts::ReconstructionType reconstruction)
{
    return SimOpts{dim,
                   0,
                   0,
                   reconstruction,
                   MHDOpts::SlopeLimiterType::None,
                   MHDOpts::RiemannSolverType::Rusanov};
}

constexpr auto Constant = MHDOpts::ReconstructionType::Constant;
constexpr auto WENOZ    = MHDOpts::ReconstructionType::WENOZ;

static_assert(mhdOpts(1, Constant).mhd_axes_consistent());
static_assert(mhdOpts(1, WENOZ).mhd_axes_consistent());


template<std::size_t dim, std::size_t interp>
struct HybridParam
{
    using Model = PHARE_Types<SimOpts{dim, interp}>::Hybrid;

    static constexpr bool is_mhd      = false;
    static constexpr auto nodalScalar = HybridQuantity::Scalar::rho;
    static constexpr auto nodalVector = HybridQuantity::Vector::V;
    static constexpr auto edgeVector  = HybridQuantity::Vector::E;
    static constexpr auto faceVector  = HybridQuantity::Vector::B;

    static constexpr auto dualScalar(std::size_t d)
    {
        return d == 0 ? HybridQuantity::Scalar::By : HybridQuantity::Scalar::Bx;
    }
};

template<std::size_t dim, MHDOpts::ReconstructionType reconstruction>
struct MHDParam
{
    using Model = PHARE_Types<mhdOpts(dim, reconstruction)>::MHD;

    static constexpr bool is_mhd      = true;
    static constexpr auto nodalScalar = MHDQuantity::Scalar::ScalarAllPrimal;
    static constexpr auto cellScalar  = MHDQuantity::Scalar::rho;
    static constexpr auto nodalVector = MHDQuantity::Vector::VecAllPrimal;
    static constexpr auto edgeVector  = MHDQuantity::Vector::E;
    static constexpr auto faceVector  = MHDQuantity::Vector::B;

    static constexpr auto dualScalar(std::size_t) { return MHDQuantity::Scalar::rho; }
};


constexpr std::array<double, 3> meshSizes{0.1, 0.2, 0.3};
constexpr std::array<std::uint32_t, 3> cellCounts{20, 14, 10};

constexpr std::array<double, 3> quadA{0.7, -1.3, 0.4};
constexpr std::array<double, 3> quadB{0.2, 0.5, -0.9};
constexpr std::array<double, 3> cubicC{0.3, -0.2, 0.15};
constexpr std::array<double, 3> waveK{1.1, 0.7, 1.9};
constexpr std::array<double, 3> wavePhase{0.0, 0.4, 1.3};

template<std::size_t dim, typename T>
std::array<T, dim> leading(std::array<T, 3> const& values)
{
    std::array<T, dim> out;
    std::copy_n(values.begin(), dim, out.begin());
    return out;
}


template<std::size_t dim>
double quadratic(std::array<double, dim> const& x)
{
    double value = 0.;
    for (std::size_t d = 0; d < dim; ++d)
        value += quadA[d] * x[d] * x[d] + quadB[d] * x[d];
    if constexpr (dim > 1)
        value += 0.6 * x[0] * x[1];
    if constexpr (dim > 2)
        value -= 0.35 * x[1] * x[2];
    return value;
}

template<std::size_t dim>
double quadraticDerivative(std::size_t d, std::array<double, dim> const& x)
{
    double value = 2. * quadA[d] * x[d] + quadB[d];
    if constexpr (dim > 1)
    {
        if (d == 0)
            value += 0.6 * x[1];
        if (d == 1)
            value += 0.6 * x[0];
    }
    if constexpr (dim > 2)
    {
        if (d == 1)
            value -= 0.35 * x[2];
        if (d == 2)
            value -= 0.35 * x[1];
    }
    return value;
}

template<std::size_t dim>
double cubic(std::array<double, dim> const& x)
{
    double value = 0.;
    for (std::size_t d = 0; d < dim; ++d)
        value += cubicC[d] * x[d] * x[d] * x[d] + quadA[d] * x[d] * x[d];
    return value;
}

template<std::size_t dim>
double cubicLaplacian(std::array<double, dim> const& x)
{
    double value = 0.;
    for (std::size_t d = 0; d < dim; ++d)
        value += 6. * cubicC[d] * x[d] + 2. * quadA[d];
    return value;
}

template<std::size_t dim>
double phase(std::array<double, dim> const& x)
{
    double value = 0.;
    for (std::size_t d = 0; d < dim; ++d)
        value += waveK[d] * x[d];
    return value;
}


template<typename Model>
struct OwnedField
{
    template<typename Layout, typename Qty>
    OwnedField(Layout const& layout, Qty qty)
        : grid{"f", layout, qty, std::numeric_limits<double>::quiet_NaN()}
        , field{"f", qty}
    {
        field.setBuffer(&grid);
    }

    OwnedField(OwnedField const&)            = delete;
    OwnedField(OwnedField&&)                 = delete;
    OwnedField& operator=(OwnedField const&) = delete;
    OwnedField& operator=(OwnedField&&)      = delete;

    Model::Grid_t grid;
    Model::Field_t field;
};

template<typename Model>
struct OwnedVecField
{
    template<typename Layout, typename Qty>
    OwnedVecField(Layout const& layout, Qty qty)
        : view{"v", qty}
        , grids{makeGrid(layout, 0), makeGrid(layout, 1), makeGrid(layout, 2)}
    {
        view.setBuffer(&grids);
    }

    OwnedVecField(OwnedVecField const&)            = delete;
    OwnedVecField(OwnedVecField&&)                 = delete;
    OwnedVecField& operator=(OwnedVecField const&) = delete;
    OwnedVecField& operator=(OwnedVecField&&)      = delete;

    Model::VecField_t view;
    std::array<typename Model::Grid_t, 3> grids;

private:
    template<typename Layout>
    auto makeGrid(Layout const& layout, std::size_t i) const
    {
        return typename Model::Grid_t{view[i].name(), layout, view[i].physicalQuantity(),
                                      std::numeric_limits<double>::quiet_NaN()};
    }
};

template<typename Grid>
double maxAbsOf(Grid const& grid)
{
    double value = 0.;
    for (auto const v : grid)
        value = std::max(value, std::abs(v));
    return value;
}

template<typename Model>
double maxAbs(OwnedField<Model> const& owned)
{
    return maxAbsOf(owned.grid);
}

template<typename Model>
double maxAbs(OwnedVecField<Model> const& owned)
{
    return std::max({maxAbsOf(owned.grids[0]), maxAbsOf(owned.grids[1]), maxAbsOf(owned.grids[2])});
}


template<typename Layout>
using Centering_t = std::array<QtyCentering, Layout::dimension>;

template<typename Layout>
using Box_t = Box<std::size_t, Layout::dimension>;

template<typename Layout>
Box_t<Layout> ghostBox(Layout const& layout, Centering_t<Layout> const& centering)
{
    Box_t<Layout> box;
    for (std::size_t d = 0; d < Layout::dimension; ++d)
    {
        box.lower[d] = layout.ghostStartIndex(centering[d], static_cast<Direction>(d));
        box.upper[d] = layout.ghostEndIndex(centering[d], static_cast<Direction>(d));
    }
    return box;
}

template<typename Layout>
Box_t<Layout> physicalBox(Layout const& layout, Centering_t<Layout> const& centering)
{
    Box_t<Layout> box;
    for (std::size_t d = 0; d < Layout::dimension; ++d)
    {
        box.lower[d] = layout.physicalStartIndex(centering[d], static_cast<Direction>(d));
        box.upper[d] = layout.physicalEndIndex(centering[d], static_cast<Direction>(d));
    }
    return box;
}

template<typename Box>
Box shrink_idx(Box box, std::size_t d)
{
    ++box.lower[d];
    --box.upper[d];
    return box;
}

template<typename Layout>
Box_t<Layout> reach(Layout const& layout, Centering_t<Layout> const& src,
                    Centering_t<Layout> const& dst, std::size_t d)
{
    auto const box = ghostBox(layout, dst);
    if (d >= Layout::dimension)
        return box;
    bool const staggered = src[d] == QtyCentering::primal and dst[d] == QtyCentering::dual;
    return staggered ? box : shrink_idx(box, d);
}

template<typename Box>
Box intersect(Box const& a, Box const& b)
{
    auto const box = a * b;
    if (!box)
        throw std::runtime_error("empty box intersection");
    return *box;
}

template<typename Layout>
std::array<double, Layout::dimension> coordinates(Layout const& layout,
                                                  Centering_t<Layout> const& centering,
                                                  MeshIndex<Layout::dimension> const& idx)
{
    auto const& dl = layout.meshSize();
    std::array<double, Layout::dimension> x;
    for (std::size_t d = 0; d < Layout::dimension; ++d)
    {
        auto const start = layout.physicalStartIndex(centering[d], static_cast<Direction>(d));
        auto const shift = centering[d] == QtyCentering::dual ? dl[d] / 2 : 0.;
        x[d] = (static_cast<double>(idx[d]) - static_cast<double>(start)) * dl[d] + shift;
    }
    return x;
}

template<typename Layout>
double tolerance(Layout const& layout, double scale)
{
    auto const& dl = layout.meshSize();
    auto const h   = *std::min_element(dl.begin(), dl.end());
    return 1e-12 * scale / (h * h);
}

template<typename Layout, typename Field, typename Fn>
void fill(Layout const& layout, Field& field, Fn const& fn)
{
    auto const centering = layout.centering(field);
    for_each_point(ghostBox(layout, centering),
                   [&](auto const& idx) { field(idx) = fn(coordinates(layout, centering, idx)); });
}

template<typename Layout, typename VecField>
void fillWave(Layout const& layout, VecField& v)
{
    for (std::size_t c = 0; c < 3; ++c)
        fill(layout, v[c], [&](auto const& x) { return std::sin(phase(x) + wavePhase[c]); });
}

template<typename Op, typename Field, typename Box>
void store(Op const& op, Field& out, Box const& box)
{
    for_each_point(box, [&](auto const& idx) { out(idx) = op(idx); });
}

template<std::size_t i, typename Vector>
auto const& component(Vector const& v)
{
    if constexpr (i == 0)
        return v.x;
    else if constexpr (i == 1)
        return v.y;
    else
        return v.z;
}

template<typename Layout, typename VecField>
void storeCurl(Layout const& layout, VecField const& v, VecField& out)
{
    auto const curlV = curl(layout, v);
    for_N<3>([&](auto ic) {
        constexpr std::size_t i = ic();
        constexpr std::size_t j = (i + 1) % 3;
        constexpr std::size_t k = (i + 2) % 3;
        auto const& op          = component<i>(curlV);
        EXPECT_EQ(op.centering(), layout.centering(out[i]));
        auto const box = intersect(reach(layout, layout.centering(v[k]), op.centering(), j),
                                   reach(layout, layout.centering(v[j]), op.centering(), k));
        store(op, out[i], box);
    });
}

template<typename Layout, typename Field, typename VecField>
void storeGrad(Layout const& layout, Field const& f, VecField& out)
{
    auto const gradF = grad(layout, f);
    for_N<3>([&](auto ic) {
        constexpr std::size_t i = ic();
        auto const& op          = component<i>(gradF);
        EXPECT_EQ(op.centering(), layout.centering(out[i]));
        store(op, out[i], reach(layout, layout.centering(f), op.centering(), i));
    });
}


template<typename Param>
struct DifferentialOperators : public ::testing::Test
{
    using Model  = Param::Model;
    using Layout = Model::GridLayout_t;

    static constexpr std::size_t dim = Layout::dimension;

    Layout layout{leading<dim>(meshSizes), leading<dim>(cellCounts), Point<double, dim>{}};
};

using DifferentialOperatorsParams
    = ::testing::Types<HybridParam<1, 1>, HybridParam<1, 3>, MHDParam<1, Constant>,
                       MHDParam<1, WENOZ>, HybridParam<2, 1>, HybridParam<2, 3>,
                       MHDParam<2, Constant>, MHDParam<2, WENOZ>, HybridParam<3, 1>,
                       HybridParam<3, 3>, MHDParam<3, Constant>, MHDParam<3, WENOZ>>;

TYPED_TEST_SUITE(DifferentialOperators, DifferentialOperatorsParams);

using LowDimParams = ::testing::Types<HybridParam<1, 1>, HybridParam<1, 3>, MHDParam<1, Constant>,
                                      MHDParam<1, WENOZ>, HybridParam<2, 1>, HybridParam<2, 3>,
                                      MHDParam<2, Constant>, MHDParam<2, WENOZ>>;

template<typename Param>
struct NonSimulatedDirections : public DifferentialOperators<Param>
{
};

TYPED_TEST_SUITE(NonSimulatedDirections, LowDimParams);


TYPED_TEST(DifferentialOperators, partialIsExactOnQuadratics)
{
    using Model               = TypeParam::Model;
    static constexpr auto dim = TestFixture::dim;
    auto const& layout        = this->layout;

    for_N<dim>([&](auto dc) {
        constexpr std::size_t d = dc();
        constexpr auto dir      = static_cast<Direction>(d);

        for (auto const qty : {TypeParam::nodalScalar, TypeParam::dualScalar(d)})
        {
            OwnedField<Model> src{layout, qty};
            fill(layout, src.field, [](auto const& x) { return quadratic(x); });
            auto const srcCentering = layout.centering(src.field);

            auto const dx = partial<dir>(layout, src.field);
            auto expected = srcCentering;
            expected[d]   = srcCentering[d] == QtyCentering::primal ? QtyCentering::dual
                                                                    : QtyCentering::primal;
            EXPECT_EQ(dx.centering(), expected);

            for_each_point(reach(layout, srcCentering, dx.centering(), d), [&](auto const& idx) {
                EXPECT_NEAR(dx(idx),
                            quadraticDerivative(d, coordinates(layout, dx.centering(), idx)),
                            1e-10);
            });
        }
    });
}


TYPED_TEST(DifferentialOperators, destinationsFollowYee)
{
    using Model               = TypeParam::Model;
    static constexpr auto dim = TestFixture::dim;
    auto const& layout        = this->layout;

    OwnedVecField<Model> e{layout, TypeParam::edgeVector};
    OwnedVecField<Model> b{layout, TypeParam::faceVector};
    OwnedField<Model> n{layout, TypeParam::nodalScalar};

    auto const curlE = curl(layout, e.view);
    auto const curlB = curl(layout, b.view);
    auto const gradN = grad(layout, n.field);
    for_N<3>([&](auto ic) {
        constexpr std::size_t i = ic();
        EXPECT_EQ(component<i>(curlE).centering(), layout.centering(b.view[i]));
        EXPECT_EQ(component<i>(curlB).centering(), layout.centering(e.view[i]));
        EXPECT_EQ(component<i>(gradN).centering(), layout.centering(e.view[i]));
    });

    if constexpr (TypeParam::is_mhd)
    {
        OwnedField<Model> cell{layout, TypeParam::cellScalar};
        auto const gradCell = grad(layout, cell.field);
        for_N<3>([&](auto ic) {
            constexpr std::size_t i = ic();
            EXPECT_EQ(component<i>(gradCell).centering(), layout.centering(b.view[i]));
        });
    }

    EXPECT_EQ(div(layout, e.view).centering(), layout.centering(n.field));
    EXPECT_EQ(div(layout, b.view).centering(), (ConstArray<QtyCentering, dim>(QtyCentering::dual)));
}


TYPED_TEST(DifferentialOperators, laplacianIsExactOnCubics)
{
    using Model        = TypeParam::Model;
    auto const& layout = this->layout;

    for (auto const qty : {TypeParam::nodalScalar, TypeParam::dualScalar(0)})
    {
        OwnedField<Model> f{layout, qty};
        fill(layout, f.field, [](auto const& x) { return cubic(x); });

        auto const centering = layout.centering(f.field);
        auto const tol       = tolerance(layout, maxAbs(f));
        for_each_point(shrink(ghostBox(layout, centering), 1), [&](auto const& idx) {
            EXPECT_NEAR(laplacian(layout, f.field, idx),
                        cubicLaplacian(coordinates(layout, centering, idx)), tol);
        });
    }
}


TYPED_TEST(DifferentialOperators, divCurlVanishes)
{
    using Model        = TypeParam::Model;
    auto const& layout = this->layout;

    auto const check = [&](auto const aQty, auto const curlQty) {
        OwnedVecField<Model> a{layout, aQty};
        fillWave(layout, a.view);
        OwnedVecField<Model> curlA{layout, curlQty};
        storeCurl(layout, a.view, curlA.view);

        auto const divA = div(layout, curlA.view);
        auto const tol  = tolerance(layout, maxAbs(a));
        for_each_point(physicalBox(layout, divA.centering()),
                       [&](auto const& idx) { EXPECT_NEAR(divA(idx), 0., tol); });
    };

    check(TypeParam::edgeVector, TypeParam::faceVector);
    check(TypeParam::faceVector, TypeParam::edgeVector);
}


TYPED_TEST(DifferentialOperators, curlGradVanishes)
{
    using Model        = TypeParam::Model;
    auto const& layout = this->layout;

    auto const check = [&](auto const phiQty, auto const gradQty) {
        OwnedField<Model> phi{layout, phiQty};
        fill(layout, phi.field, [](auto const& x) { return std::sin(phase(x)); });
        OwnedVecField<Model> gradPhi{layout, gradQty};
        storeGrad(layout, phi.field, gradPhi.view);

        auto const tol      = tolerance(layout, maxAbs(phi));
        auto const curlGrad = curl(layout, gradPhi.view);
        for_N<3>([&](auto ic) {
            constexpr std::size_t i = ic();
            auto const& op          = component<i>(curlGrad);
            for_each_point(physicalBox(layout, op.centering()),
                           [&](auto const& idx) { EXPECT_NEAR(op(idx), 0., tol); });
        });
    };

    check(TypeParam::nodalScalar, TypeParam::edgeVector);
    if constexpr (TypeParam::is_mhd)
        check(TypeParam::cellScalar, TypeParam::faceVector);
}


TYPED_TEST(DifferentialOperators, divGradEqualsLaplacian)
{
    using Model        = TypeParam::Model;
    auto const& layout = this->layout;

    auto const check = [&](auto const phiQty, auto const gradQty) {
        OwnedField<Model> phi{layout, phiQty};
        fill(layout, phi.field, [](auto const& x) { return std::sin(phase(x)); });
        OwnedVecField<Model> gradPhi{layout, gradQty};
        storeGrad(layout, phi.field, gradPhi.view);

        auto const divGrad = div(layout, gradPhi.view);
        EXPECT_EQ(divGrad.centering(), layout.centering(phi.field));
        auto const tol = tolerance(layout, maxAbs(phi));
        for_each_point(physicalBox(layout, layout.centering(phi.field)), [&](auto const& idx) {
            EXPECT_NEAR(divGrad(idx), laplacian(layout, phi.field, idx), tol);
        });
    };

    check(TypeParam::nodalScalar, TypeParam::edgeVector);
    if constexpr (TypeParam::is_mhd)
        check(TypeParam::cellScalar, TypeParam::faceVector);
}


TYPED_TEST(NonSimulatedDirections, nonSimulatedDirectionsAreZero)
{
    using Model               = TypeParam::Model;
    static constexpr auto dim = TestFixture::dim;
    auto const& layout        = this->layout;

    OwnedField<Model> f{layout, TypeParam::nodalScalar};
    fill(layout, f.field, [](auto const& x) { return quadratic(x); });
    auto const nodal = layout.centering(f.field);

    auto const dz = partial<Direction::Z>(layout, f.field);
    EXPECT_EQ(dz.centering(), nodal);
    for_each_point(ghostBox(layout, nodal), [&](auto const& idx) { EXPECT_EQ(dz(idx), 0.); });

    if constexpr (dim == 1)
    {
        auto const dy = partial<Direction::Y>(layout, f.field);
        EXPECT_EQ(dy.centering(), nodal);
        for_each_point(ghostBox(layout, nodal), [&](auto const& idx) { EXPECT_EQ(dy(idx), 0.); });

        OwnedVecField<Model> a{layout, TypeParam::edgeVector};
        fillWave(layout, a.view);

        auto const curlA = curl(layout, a.view);
        auto const dxAz  = partial<Direction::X>(layout, a.view[2]);
        auto const dxAy  = partial<Direction::X>(layout, a.view[1]);

        for_each_point(ghostBox(layout, curlA.x.centering()), [&](auto const& idx) {
            EXPECT_EQ(curlA.x(idx), 0.);
            EXPECT_FALSE(std::signbit(curlA.x(idx)));
        });
        for_each_point(ghostBox(layout, curlA.y.centering()), [&](auto const& idx) {
            EXPECT_EQ(curlA.y(idx), -dxAz(idx));
            EXPECT_EQ(std::signbit(curlA.y(idx)), std::signbit(-dxAz(idx)));
        });
        for_each_point(ghostBox(layout, curlA.z.centering()), [&](auto const& idx) {
            EXPECT_EQ(curlA.z(idx), dxAy(idx));
            EXPECT_EQ(std::signbit(curlA.z(idx)), std::signbit(dxAy(idx)));
        });
    }
}


TYPED_TEST(DifferentialOperators, collocatedVectorThrows)
{
    using Model        = TypeParam::Model;
    auto const& layout = this->layout;

    OwnedVecField<Model> v{layout, TypeParam::nodalVector};
    EXPECT_THROW(curl(layout, v.view), std::runtime_error);
    EXPECT_THROW(div(layout, v.view), std::runtime_error);
}


TYPED_TEST(DifferentialOperators, operandWithoutDataThrows)
{
    using Model        = TypeParam::Model;
    auto const& layout = this->layout;

    typename Model::Field_t f{"f", TypeParam::nodalScalar};
    EXPECT_THROW(partial<Direction::X>(layout, f), std::runtime_error);
    EXPECT_THROW(grad(layout, f), std::runtime_error);
}


TYPED_TEST(DifferentialOperators, componentArrayAndVecFieldAgree)
{
    using Model        = TypeParam::Model;
    auto const& layout = this->layout;

    OwnedVecField<Model> e{layout, TypeParam::edgeVector};
    fillWave(layout, e.view);

    std::array const components{e.view[0], e.view[1], e.view[2]};
    auto const fromVecField   = curl(layout, e.view);
    auto const fromComponents = curl(layout, components);
    for_N<3>([&](auto ic) {
        constexpr std::size_t i = ic();
        auto const& expected    = component<i>(fromVecField);
        auto const& actual      = component<i>(fromComponents);
        for_each_point(ghostBox(layout, expected.centering()), [&](auto const& idx) {
            EXPECT_EQ(std::bit_cast<std::uint64_t>(expected(idx)),
                      std::bit_cast<std::uint64_t>(actual(idx)));
        });
    });
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
