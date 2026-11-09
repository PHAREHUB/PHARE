#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/numerics/constrained_transport/non_ideal_electric_field.hpp"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace PHARE::core
{

template<std::size_t dim>
struct MHDParam
{
    static constexpr SimOpts opts{dim,
                                  0,
                                  0,
                                  MHDOpts::ReconstructionType::Constant,
                                  MHDOpts::SlopeLimiterType::None,
                                  MHDOpts::RiemannSolverType::Rusanov};
    using Model = PHARE_Types<opts>::MHD;
};

constexpr std::array<double, 3> meshSizes{0.1, 0.2, 0.3};
constexpr std::array<std::uint32_t, 3> cellCounts{20, 14, 10};

template<std::size_t dim, typename T>
std::array<T, dim> leading(std::array<T, 3> const& values)
{
    std::array<T, dim> out;
    std::copy_n(values.begin(), dim, out.begin());
    return out;
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
    OwnedField& operator=(OwnedField const&) = delete;

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
    OwnedVecField& operator=(OwnedVecField const&) = delete;

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

template<typename Layout>
std::array<double, 3> coordinates(Layout const& layout, Centering_t<Layout> const& centering,
                                  MeshIndex<Layout::dimension> const& idx)
{
    auto const& dl = layout.meshSize();
    std::array<double, 3> x{0., 0., 0.};
    for (std::size_t d = 0; d < Layout::dimension; ++d)
    {
        auto const start = layout.physicalStartIndex(centering[d], static_cast<Direction>(d));
        auto const shift = centering[d] == QtyCentering::dual ? dl[d] / 2 : 0.;
        x[d] = (static_cast<double>(idx[d]) - static_cast<double>(start)) * dl[d] + shift;
    }
    return x;
}

template<typename Layout, typename Field, typename Fn>
void fill(Layout const& layout, Field& field, Fn const& fn)
{
    auto const centering = layout.centering(field);
    for_each_point(ghostBox(layout, centering),
                   [&](auto const& idx) { field(idx) = fn(coordinates(layout, centering, idx)); });
}

template<typename VecField>
struct NonIdealState
{
    VecField& e;
    auto& E() { return e; }
    auto const& E() const { return e; }
};

template<typename VecField, typename Field>
struct MHDView
{
    VecField& J;
    VecField& B;
    Field& rho;
};


template<typename Param>
struct NonIdealElectricFieldTest : public ::testing::Test
{
    using Model                      = Param::Model;
    using Layout                     = Model::GridLayout_t;
    static constexpr std::size_t dim = Layout::dimension;

    Layout layout{leading<dim>(meshSizes), leading<dim>(cellCounts), Point<double, dim>{}};

    OwnedVecField<Model> J{layout, MHDQuantity::Vector::J};
    OwnedVecField<Model> B{layout, MHDQuantity::Vector::B};
    OwnedVecField<Model> E{layout, MHDQuantity::Vector::E};
    OwnedField<Model> rho{layout, MHDQuantity::Scalar::rho};

    template<bool Hall>
    void compute(OhmInfo const& info)
    {
        NonIdealState<typename Model::VecField_t> state{E.view};
        MHDView<typename Model::VecField_t, typename Model::Field_t> mhd{J.view, B.view, rho.field};
        NonIdealElectricField<Layout, Hall>{info, layout}(state, mhd);
    }
};

using NonIdealParams = ::testing::Types<MHDParam<1>, MHDParam<2>, MHDParam<3>>;

TYPED_TEST_SUITE(NonIdealElectricFieldTest, NonIdealParams);


double affine(std::array<double, 3> const& x, double a, double bx, double by, double bz)
{
    return a + bx * x[0] + by * x[1] + bz * x[2];
}


TYPED_TEST(NonIdealElectricFieldTest, resistiveOnlyIsEtaJ)
{
    auto const& layout = this->layout;
    double const eta   = 0.3;

    for (std::size_t c = 0; c < 3; ++c)
    {
        fill(layout, this->J.view[c],
             [&](auto const& x) { return affine(x, 0.4 + c, 0.7, -0.2 * c, 0.3); });
        fill(layout, this->B.view[c],
             [&](auto const& x) { return affine(x, 1. + c, 0.1, 0., 0.); });
    }
    fill(layout, this->rho.field, [&](auto const& x) { return affine(x, 2., 0.1, 0.2, -0.1); });

    this->template compute<false>(OhmInfo{eta, 0.0, HyperMode::constant});

    for (std::size_t c = 0; c < 3; ++c)
    {
        auto const& Ec = this->E.view[c];
        auto const& Jc = this->J.view[c];
        for_each_point(physicalBox(layout, layout.centering(Ec)),
                       [&](auto const& idx) { EXPECT_DOUBLE_EQ(Ec(idx), eta * Jc(idx)); });
    }
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
