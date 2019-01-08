#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "data/field/field.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/grid/gridlayoutdefs.h"
#include "data/vecfield/vecfield.h"
#include "numerics/faraday/faraday.h"
#include "utilities/box/box.h"
#include "utilities/index/index.h"
#include "utilities/point/point.h"

using namespace PHARE;

struct FieldMock
{
    double data;
    double& operator()(uint32 i) { return data; }
    double const& operator()(uint32 i) const { return data; }
    QtyCentering physicalQuantity() { return QtyCentering::dual; }
};

struct VecFieldMock
{
    FieldMock fm;
    FieldMock& getComponent(Component comp) { return fm; }
    FieldMock const& getComponent(Component comp) const { return fm; }
};


struct GridLayoutMock1D
{
    static const int dimension = 1;
    double deriv(FieldMock const& f, MeshIndex<1> mi, DirectionTag<Direction::X>) {}
    int physicalStartIndex(FieldMock&, Direction dir) { return 0; }
    int physicalEndIndex(FieldMock&, Direction dir) { return 0; }
};

struct GridLayoutMock2D
{
    static const int dimension = 2;
    double deriv(FieldMock const& f, MeshIndex<2> mi, DirectionTag<Direction::X>) { return 0; }
    double deriv(FieldMock const& f, MeshIndex<2> mi, DirectionTag<Direction::Y>) { return 0; }
    int physicalStartIndex(FieldMock&, Direction dir) { return 0; }
    int physicalEndIndex(FieldMock&, Direction dir) { return 0; }
};

struct GridLayoutMock3D
{
    static const int dimension = 3;
    double deriv(FieldMock const& f, MeshIndex<3> mi, DirectionTag<Direction::X>) { return 0; }
    double deriv(FieldMock const& f, MeshIndex<3> mi, DirectionTag<Direction::Y>) { return 0; }
    double deriv(FieldMock const& f, MeshIndex<3> mi, DirectionTag<Direction::Z>) { return 0; }
    int physicalStartIndex(FieldMock&, Direction dir) { return 0; }
    int physicalEndIndex(FieldMock&, Direction dir) { return 0; }
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
    VecFieldMock B, E, Bnew;

    Faraday<GridLayoutMock1D> faraday1d;
    auto layout1d = std::make_unique<GridLayoutMock1D>();
    EXPECT_ANY_THROW(faraday1d(B, E, Bnew));
    faraday1d.setLayout(layout1d.get());

    Faraday<GridLayoutMock2D> faraday2d;
    auto layout2d = std::make_unique<GridLayoutMock2D>();
    // EXPECT_ANY_THROW(faraday2d(B, E, Bnew));
    faraday2d.setLayout(layout2d.get());

    Faraday<GridLayoutMock3D> faraday3d;
    auto layout3d = std::make_unique<GridLayoutMock3D>();
    // EXPECT_ANY_THROW(faraday3d(B, E, Bnew));
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
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Bx;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> By;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Bz;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Ex;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Ey;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Ez;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Bxnew;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Bynew;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Bznew;
    VecField<NdArrayVector1D<>, HybridQuantity> B;
    VecField<NdArrayVector1D<>, HybridQuantity> E;
    VecField<NdArrayVector1D<>, HybridQuantity> Bnew;
    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday1DTest()
        : layout{{{0.1}}, {{50}}, Point{0.}, Box{Point{0}, Point{49}}}
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
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> Bx;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> By;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> Bz;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> Ex;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> Ey;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> Ez;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> Bxnew;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> Bynew;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> Bznew;
    VecField<NdArrayVector2D<>, HybridQuantity> B;
    VecField<NdArrayVector2D<>, HybridQuantity> E;
    VecField<NdArrayVector2D<>, HybridQuantity> Bnew;
    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday2DTest()
        : layout{{{0.1, 0.2}}, {{50, 30}}, Point{0., 0.}, Box{Point{0, 0}, Point{49, 39}}}
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
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> Bx;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> By;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> Bz;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> Ex;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> Ey;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> Ez;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> Bxnew;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> Bynew;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> Bznew;
    VecField<NdArrayVector3D<>, HybridQuantity> B;
    VecField<NdArrayVector3D<>, HybridQuantity> E;
    VecField<NdArrayVector3D<>, HybridQuantity> Bnew;
    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday3DTest()
        : layout{{{0.1, 0.2, 0.3}},
                 {{50, 30, 40}},
                 Point{0., 0., 0.},
                 Box{Point{0, 0, 0}, Point{49, 29, 39}}}
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
    faraday(B, E, Bnew);

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
    faraday(B, E, Bnew);

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
    faraday(B, E, Bnew);

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
