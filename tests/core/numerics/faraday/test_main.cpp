#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "core/data/field/field.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield.h"
#include "core/numerics/faraday/faraday.h"
#include "core/utilities/box/box.h"
#include "core/utilities/index/index.h"
#include "core/utilities/point/point.h"

using namespace PHARE::core;

template<std::size_t dim>
struct FieldMock
{
    static auto constexpr dimension = dim;
    double data;
    double& operator()([[maybe_unused]] std::uint32_t i) { return data; }
    double const& operator()([[maybe_unused]] std::uint32_t i) const { return data; }
    double& operator()([[maybe_unused]] std::uint32_t i, [[maybe_unused]] std::uint32_t j)
    {
        return data;
    }
    double const& operator()([[maybe_unused]] std::uint32_t i,
                             [[maybe_unused]] std::uint32_t j) const
    {
        return data;
    }
    double& operator()([[maybe_unused]] std::uint32_t i, [[maybe_unused]] std::uint32_t j,
                       [[maybe_unused]] std::uint32_t k)
    {
        return data;
    }
    double const& operator()([[maybe_unused]] std::uint32_t i, [[maybe_unused]] std::uint32_t j,
                             [[maybe_unused]] std::uint32_t k) const
    {
        return data;
    }
    QtyCentering physicalQuantity() { return QtyCentering::dual; }
    std::string name() const { return "FieldMock"; }
};

template<typename Field>
struct VecFieldMock
{
    using field_type                = Field;
    static auto constexpr dimension = Field::dimension;
    Field fm;
    Field& getComponent([[maybe_unused]] Component comp) { return fm; }
    Field const& getComponent([[maybe_unused]] Component comp) const { return fm; }
    bool isUsable() const { return true; }
};


struct GridLayoutMock1D
{
    static const auto dimension = 1u;
    double deriv([[maybe_unused]] FieldMock<1> const& f, [[maybe_unused]] MeshIndex<1u> mi,
                 [[maybe_unused]] DirectionTag<Direction::X>)
    {
        return 0;
    }
    std::size_t physicalStartIndex([[maybe_unused]] FieldMock<1>&, [[maybe_unused]] Direction dir)
    {
        return 0;
    }
    std::size_t physicalEndIndex([[maybe_unused]] FieldMock<1>&, [[maybe_unused]] Direction dir)
    {
        return 0;
    }
};

struct GridLayoutMock2D
{
    static const auto dimension = 2u;
    double deriv([[maybe_unused]] FieldMock<dimension> const& f, [[maybe_unused]] MeshIndex<2u> mi,
                 [[maybe_unused]] DirectionTag<Direction::X>)
    {
        return 0;
    }
    double deriv([[maybe_unused]] FieldMock<dimension> const& f, [[maybe_unused]] MeshIndex<2u> mi,
                 [[maybe_unused]] DirectionTag<Direction::Y>)
    {
        return 0;
    }
    std::size_t physicalStartIndex([[maybe_unused]] FieldMock<dimension>&,
                                   [[maybe_unused]] Direction dir)
    {
        return 0;
    }
    std::size_t physicalEndIndex([[maybe_unused]] FieldMock<dimension>&,
                                 [[maybe_unused]] Direction dir)
    {
        return 0;
    }
};

struct GridLayoutMock3D
{
    static const auto dimension = 3u;
    double deriv([[maybe_unused]] FieldMock<dimension> const& f, [[maybe_unused]] MeshIndex<3u> mi,
                 [[maybe_unused]] DirectionTag<Direction::X>)
    {
        return 0;
    }
    double deriv([[maybe_unused]] FieldMock<dimension> const& f, [[maybe_unused]] MeshIndex<3u> mi,
                 [[maybe_unused]] DirectionTag<Direction::Y>)
    {
        return 0;
    }
    double deriv([[maybe_unused]] FieldMock<dimension> const& f, [[maybe_unused]] MeshIndex<3u> mi,
                 [[maybe_unused]] DirectionTag<Direction::Z>)
    {
        return 0;
    }
    std::size_t physicalStartIndex([[maybe_unused]] FieldMock<dimension>&,
                                   [[maybe_unused]] Direction dir)
    {
        return 0;
    }
    std::size_t physicalEndIndex([[maybe_unused]] FieldMock<dimension>&,
                                 [[maybe_unused]] Direction dir)
    {
        return 0;
    }
};



TEST(Faraday, canBe1D)
{
    Faraday<GridLayoutMock1D> faraday;
}


TEST(Faraday, canBe2D)
{
    Faraday<GridLayoutMock2D> faraday;
}


TEST(Faraday, canBe3D)
{
    Faraday<GridLayoutMock3D> faraday;
}


TEST(Faraday, shouldBeGivenAGridLayoutPointerToBeOperational)
{
    VecFieldMock<FieldMock<1>> B_1, E_1, Bnew_1;

    Faraday<GridLayoutMock1D> faraday1d;
    auto layout1d = std::make_unique<GridLayoutMock1D>();
    EXPECT_ANY_THROW(faraday1d(B_1, E_1, Bnew_1, 1.));
    faraday1d.setLayout(layout1d.get());

    VecFieldMock<FieldMock<2>> B_2, E_2, Bnew_2;

    Faraday<GridLayoutMock2D> faraday2d;
    auto layout2d = std::make_unique<GridLayoutMock2D>();
    EXPECT_ANY_THROW(faraday2d(B_2, E_2, Bnew_2, 1.));
    faraday2d.setLayout(layout2d.get());

    VecFieldMock<FieldMock<3>> B_3, E_3, Bnew_3;

    Faraday<GridLayoutMock3D> faraday3d;
    auto layout3d = std::make_unique<GridLayoutMock3D>();
    EXPECT_ANY_THROW(faraday3d(B_3, E_3, Bnew_3, 1.));
    faraday3d.setLayout(layout3d.get());
}




std::vector<double> read(std::string filename)
{
    std::ifstream readFile(filename);
    assert(readFile.is_open());
    std::vector<double> x;

    std::copy(std::istream_iterator<double>(readFile), std::istream_iterator<double>(),
              std::back_inserter(x));
    return x;
}



class Faraday1DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<1, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr auto interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Ex;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Ey;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Ez;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bxnew;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bynew;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bznew;
    VecField<NdArrayVector<1>, HybridQuantity> B;
    VecField<NdArrayVector<1>, HybridQuantity> E;
    VecField<NdArrayVector<1>, HybridQuantity> Bnew;
    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday1DTest()
        : layout{{{0.1}}, {{50}}, Point{0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"Ex", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Ey{"Ey", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Ez{"Ez", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , Bxnew{"Bxnew", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , Bynew{"Bynew", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bznew{"Bznew", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , B{"B", HybridQuantity::Vector::B}
        , E{"E", HybridQuantity::Vector::E}
        , Bnew{"Bnew", HybridQuantity::Vector::B}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        E.setBuffer("E_x", &Ex);
        E.setBuffer("E_y", &Ey);
        E.setBuffer("E_z", &Ez);
        Bnew.setBuffer("Bnew_x", &Bxnew);
        Bnew.setBuffer("Bnew_y", &Bynew);
        Bnew.setBuffer("Bnew_z", &Bznew);
    }
};




class Faraday2DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<2, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr auto interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Ex;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Ey;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Ez;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bxnew;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bynew;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bznew;
    VecField<NdArrayVector<2>, HybridQuantity> B;
    VecField<NdArrayVector<2>, HybridQuantity> E;
    VecField<NdArrayVector<2>, HybridQuantity> Bnew;
    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday2DTest()
        : layout{{{0.1, 0.2}}, {{50, 30}}, Point{0., 0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"Ex", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Ey{"Ey", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Ez{"Ez", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , Bxnew{"Bxnew", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , Bynew{"Bynew", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bznew{"Bznew", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , B{"B", HybridQuantity::Vector::B}
        , E{"E", HybridQuantity::Vector::E}
        , Bnew{"Bnew", HybridQuantity::Vector::B}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        E.setBuffer("E_x", &Ex);
        E.setBuffer("E_y", &Ey);
        E.setBuffer("E_z", &Ez);
        Bnew.setBuffer("Bnew_x", &Bxnew);
        Bnew.setBuffer("Bnew_y", &Bynew);
        Bnew.setBuffer("Bnew_z", &Bznew);
    }
};




class Faraday3DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<3, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr auto interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Ex;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Ey;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Ez;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bxnew;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bynew;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bznew;
    VecField<NdArrayVector<3>, HybridQuantity> B;
    VecField<NdArrayVector<3>, HybridQuantity> E;
    VecField<NdArrayVector<3>, HybridQuantity> Bnew;
    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday3DTest()
        : layout{{{0.1, 0.2, 0.3}}, {{50, 30, 40}}, Point{0., 0., 0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"Ex", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Ey{"Ey", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Ez{"Ez", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , Bxnew{"Bxnew", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , Bynew{"Bynew", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bznew{"Bznew", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , B{"B", HybridQuantity::Vector::B}
        , E{"E", HybridQuantity::Vector::E}
        , Bnew{"Bnew", HybridQuantity::Vector::B}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        E.setBuffer("E_x", &Ex);
        E.setBuffer("E_y", &Ey);
        E.setBuffer("E_z", &Ez);
        Bnew.setBuffer("Bnew_x", &Bxnew);
        Bnew.setBuffer("Bnew_y", &Bynew);
        Bnew.setBuffer("Bnew_z", &Bznew);
    }
};




TEST_F(Faraday1DTest, Faraday1DCalculatedOk)
{
    auto filename_dbydt = std::string{"dbydt_yee_1D_order1.txt"};
    auto filename_dbzdt = std::string{"dbzdt_yee_1D_order1.txt"};
    auto expected_dbydt = read(filename_dbydt);
    auto expected_dbzdt = read(filename_dbzdt);

    auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        auto point = this->layout.fieldNodeCoordinates(Ey, Point{0.}, ix);

        Ey(ix) = std::cos(2 * M_PI / 5. * point[0]);
        Ez(ix) = std::sin(2 * M_PI / 5. * point[0]);
    }

    auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        auto point = this->layout.fieldNodeCoordinates(By, Point{0.}, ix);

        By(ix) = std::tanh(point[0] - 5. / 2.);
        Bz(ix) = std::tanh(point[0] - 5. / 2.);
    }

    faraday.setLayout(&layout);
    faraday(B, E, Bnew, 1.);

    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        EXPECT_THAT(Bynew(ix), ::testing::DoubleNear((expected_dbydt[ix]), 1e-12));
        EXPECT_THAT(Bznew(ix), ::testing::DoubleNear((expected_dbzdt[ix]), 1e-12));
    }
}




TEST_F(Faraday2DTest, Faraday2DCalculatedOk)
{
    auto filename_dbxdt = std::string{"dbxdt_yee_2D_order1.txt"};
    auto filename_dbydt = std::string{"dbydt_yee_2D_order1.txt"};
    auto filename_dbzdt = std::string{"dbzdt_yee_2D_order1.txt"};
    auto expected_dbxdt = read(filename_dbxdt);
    auto expected_dbydt = read(filename_dbydt);
    auto expected_dbzdt = read(filename_dbzdt);

    auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    auto gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
    auto gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Ex, Point{0., 0.}, ix, iy);

            Ex(ix, iy) = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Ey, Point{0., 0.}, ix, iy);

            Ey(ix, iy) = std::cos(2 * M_PI / 5. * point[0]) * std::tanh(2 * M_PI / 6. * point[1]);
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Ez, Point{0., 0.}, ix, iy);

            Ez(ix, iy) = std::sin(2 * M_PI / 5. * point[0]) * std::tanh(2 * M_PI / 6. * point[1]);
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Bx, Point{0., 0.}, ix, iy);

            Bx(ix, iy) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.);
        }
    }

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(By, Point{0., 0.}, ix, iy);

            By(ix, iy) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.);
        }
    }

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Bz, Point{0., 0.}, ix, iy);

            Bz(ix, iy) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.);
        }
    }

    faraday.setLayout(&layout);
    faraday(B, E, Bnew, 1.);

    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

    auto nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Bx);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(Bxnew(ix, iy), ::testing::DoubleNear((expected_dbxdt[index_]), 1e-12));
        }
    }

    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::By);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(Bynew(ix, iy), ::testing::DoubleNear((expected_dbydt[index_]), 1e-12));
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Bz);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(Bznew(ix, iy), ::testing::DoubleNear((expected_dbzdt[index_]), 1e-12));
        }
    }
}



TEST_F(Faraday3DTest, Faraday3DCalculatedOk)
{
    auto filename_dbxdt = std::string{"dbxdt_yee_3D_order1.txt"};
    auto filename_dbydt = std::string{"dbydt_yee_3D_order1.txt"};
    auto filename_dbzdt = std::string{"dbzdt_yee_3D_order1.txt"};
    auto expected_dbxdt = read(filename_dbxdt);
    auto expected_dbydt = read(filename_dbydt);
    auto expected_dbzdt = read(filename_dbzdt);

    auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    auto gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
    auto gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);

    auto gsi_p_Z = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Z);
    auto gei_p_Z = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Z);
    auto gsi_d_Z = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Z);
    auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Ex, Point{0., 0., 0.}, ix, iy, iz);

                Ex(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])
                                 * std::cos(2 * M_PI / 6. * point[1])
                                 * std::tanh(2 * M_PI / 12. * point[2]);
            }
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Ey, Point{0., 0., 0.}, ix, iy, iz);

                Ey(ix, iy, iz) = std::tanh(2 * M_PI / 5. * point[0])
                                 * std::sin(2 * M_PI / 6. * point[1])
                                 * std::cos(2 * M_PI / 12. * point[2]);
            }
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Ez, Point{0., 0., 0.}, ix, iy, iz);

                Ez(ix, iy, iz) = std::cos(2 * M_PI / 5. * point[0])
                                 * std::tanh(2 * M_PI / 6. * point[1])
                                 * std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Bx, Point{0., 0., 0.}, ix, iy, iz);

                Bx(ix, iy, iz) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.)
                                 * std::tanh(point[2] - 12. / 2.);
            }
        }
    }

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(By, Point{0., 0., 0.}, ix, iy, iz);

                By(ix, iy, iz) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.)
                                 * std::tanh(point[2] - 12. / 2.);
            }
        }
    }

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Bz, Point{0., 0., 0.}, ix, iy, iz);

                Bz(ix, iy, iz) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.)
                                 * std::tanh(point[2] - 12. / 2.);
            }
        }
    }

    faraday.setLayout(&layout);
    faraday(B, E, Bnew, 1.);

    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);

    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

    auto psi_p_Z = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Z);
    auto pei_p_Z = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Z);
    auto psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    auto pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    auto nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Bx);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            for (auto iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(Bxnew(ix, iy, iz),
                            ::testing::DoubleNear((expected_dbxdt[index_]), 1e-12));
            }
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::By);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (auto iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(Bynew(ix, iy, iz),
                            ::testing::DoubleNear((expected_dbydt[index_]), 1e-12));
            }
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Bz);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            for (auto iz = psi_p_Z; iz <= pei_p_Z; ++iz)
            {
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(Bznew(ix, iy, iz),
                            ::testing::DoubleNear((expected_dbzdt[index_]), 1e-12));
            }
        }
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
