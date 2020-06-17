#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "core/data/field/field.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield.h"
#include "core/numerics/ohm/ohm.h"
#include "core/utilities/index/index.h"

using namespace PHARE::core;

template<std::size_t dim>
struct FieldMock
{
    static auto constexpr dimension = dim;
    double data;
    double& operator()([[maybe_unused]] uint32 i) { return data; }
    double const& operator()([[maybe_unused]] uint32 i) const { return data; }
    QtyCentering physicalQuantity() { return QtyCentering::dual; }
};

template<typename Field>
struct VecFieldMock
{
    using field_type                = Field;
    static auto constexpr dimension = Field::dimension;
    Field fm;
    Field& getComponent([[maybe_unused]] Component comp) { return fm; }
    Field const& getComponent([[maybe_unused]] Component comp) const { return fm; }
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



TEST(Ohm, canBe1D)
{
    Ohm<GridLayoutMock1D> ohm;
}


TEST(Ohm, canBe2D)
{
    Ohm<GridLayoutMock2D> ohm;
}


TEST(Ohm, canBe3D)
{
    Ohm<GridLayoutMock3D> ohm;
}




TEST(Ohm, shouldBeGivenAGridLayoutPointerToBeOperational)
{
    [[maybe_unused]] FieldMock<1> n, Pe;
    [[maybe_unused]] VecFieldMock<FieldMock<1>> Ve, B, J, Enew;

    Ohm<GridLayoutMock1D> ohm1d;
    auto layout1d = std::make_unique<GridLayoutMock1D>();
    // EXPECT_ANY_THROW(ohm1d(n, Ve, Pe, B, J, Enew));// TODO issue #3392
    ohm1d.setLayout(layout1d.get());

    Ohm<GridLayoutMock2D> ohm2d;
    auto layout2d = std::make_unique<GridLayoutMock2D>();
    // EXPECT_ANY_THROW(ohm2d(n, Ve, Pe, B, J, Enew));// TODO issue #3392
    ohm2d.setLayout(layout2d.get());

    Ohm<GridLayoutMock3D> ohm3d;
    auto layout3d = std::make_unique<GridLayoutMock3D>();
    // EXPECT_ANY_THROW(ohm3d(n, Ve, Pe, B, J, Enew));// TODO issue #3392
    ohm3d.setLayout(layout3d.get());
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



class Ohm1DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<1, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> n;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Vx;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Vy;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Vz;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> P;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jz;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Exnew;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Eynew;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Eznew;
    VecField<NdArrayVector<1>, HybridQuantity> V;
    VecField<NdArrayVector<1>, HybridQuantity> B;
    VecField<NdArrayVector<1>, HybridQuantity> J;
    VecField<NdArrayVector<1>, HybridQuantity> Enew;
    Ohm<GridLayout<GridLayoutImpl>> ohm;

public:
    Ohm1DTest()
        : layout{{{0.1}}, {{50}}, Point<double, 1>{0.}}
        , n{"n", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho)}
        , Vx{"Vx", HybridQuantity::Scalar::Vx, layout.allocSize(HybridQuantity::Scalar::Vx)}
        , Vy{"Vy", HybridQuantity::Scalar::Vy, layout.allocSize(HybridQuantity::Scalar::Vy)}
        , Vz{"Vz", HybridQuantity::Scalar::Vz, layout.allocSize(HybridQuantity::Scalar::Vz)}
        , P{"P", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , Exnew{"Exnew", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Eynew{"Eynew", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Eznew{"Eznew", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
        , V{"V", HybridQuantity::Vector::V}
        , B{"B", HybridQuantity::Vector::B}
        , J{"J", HybridQuantity::Vector::J}
        , Enew{"Enew", HybridQuantity::Vector::E}
    {
        V.setBuffer("V_x", &Vx);
        V.setBuffer("V_y", &Vy);
        V.setBuffer("V_z", &Vz);
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        J.setBuffer("J_x", &Jx);
        J.setBuffer("J_y", &Jy);
        J.setBuffer("J_z", &Jz);
        Enew.setBuffer("Enew_x", &Exnew);
        Enew.setBuffer("Enew_y", &Eynew);
        Enew.setBuffer("Enew_z", &Eznew);

        auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
        auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);

        auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
        auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

        for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
        {
            auto point = this->layout.fieldNodeCoordinates(Vx, Point<double, 1>{0.}, ix);

            Vx(ix) = std::cosh(0.2 * point[0]);
            Vy(ix) = std::cosh(0.4 * point[0]);
            Vz(ix) = std::cosh(0.3 * point[0]);
            Bx(ix) = std::tanh(0.2 * point[0]);
            n(ix)  = std::sinh(0.5 * point[0]);
            P(ix)  = std::sinh(0.8 * point[0]);
            Jy(ix) = std::sinh(0.3 * point[0]);
            Jz(ix) = std::sinh(0.2 * point[0]);
        }

        for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
        {
            auto point = this->layout.fieldNodeCoordinates(By, Point<double, 1>{0.}, ix);

            By(ix) = std::tanh(3.0 * point[0]);
            Bz(ix) = std::tanh(1.2 * point[0]);
            Jx(ix) = std::sinh(0.1 * point[0]);
        }
    }
};



class Ohm2DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<2, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> n;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Vx;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Vy;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Vz;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> P;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Jz;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Exnew;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Eynew;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Eznew;
    VecField<NdArrayVector<2>, HybridQuantity> V;
    VecField<NdArrayVector<2>, HybridQuantity> B;
    VecField<NdArrayVector<2>, HybridQuantity> J;
    VecField<NdArrayVector<2>, HybridQuantity> Enew;
    Ohm<GridLayout<GridLayoutImpl>> ohm;

public:
    Ohm2DTest()
        : layout{{{0.1, 0.2}}, {{50, 30}}, Point<double, 2>{0., 0.}}
        , n{"n", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho)}
        , Vx{"Vx", HybridQuantity::Scalar::Vx, layout.allocSize(HybridQuantity::Scalar::Vx)}
        , Vy{"Vy", HybridQuantity::Scalar::Vy, layout.allocSize(HybridQuantity::Scalar::Vy)}
        , Vz{"Vz", HybridQuantity::Scalar::Vz, layout.allocSize(HybridQuantity::Scalar::Vz)}
        , P{"P", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , Exnew{"Exnew", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Eynew{"Eynew", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Eznew{"Eznew", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
        , V{"V", HybridQuantity::Vector::V}
        , B{"B", HybridQuantity::Vector::B}
        , J{"J", HybridQuantity::Vector::J}
        , Enew{"Enew", HybridQuantity::Vector::E}
    {
        V.setBuffer("V_x", &Vx);
        V.setBuffer("V_y", &Vy);
        V.setBuffer("V_z", &Vz);
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        J.setBuffer("J_x", &Jx);
        J.setBuffer("J_y", &Jy);
        J.setBuffer("J_z", &Jz);
        Enew.setBuffer("Enew_x", &Exnew);
        Enew.setBuffer("Enew_y", &Eynew);
        Enew.setBuffer("Enew_z", &Eznew);

        auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
        auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
        auto gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
        auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);

        auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
        auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
        auto gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
        auto gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);

        for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
        {
            for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
            {
                auto point = this->layout.fieldNodeCoordinates(
                    Vx, Point<double, 2>{0., 0.}, ix,
                    iy); // only 1 point as all quantities are primal,primal

                Vx(ix, iy) = std::cosh(0.2 * point[0]) * std::sinh(0.2 * point[1]);
                Vy(ix, iy) = std::cosh(0.4 * point[0]) * std::sinh(0.4 * point[1]);
                Vz(ix, iy) = std::cosh(0.3 * point[0]) * std::sinh(0.3 * point[1]);
                n(ix, iy)  = std::exp(-0.1 * point[0]) * std::exp(-0.1 * point[1]);
                P(ix, iy)  = std::exp(-0.2 * point[0]) * std::exp(-0.2 * point[1]);
                Jz(ix, iy) = std::sinh(0.2 * point[0]) * std::sinh(0.2 * point[1]);
            }

            for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
            {
                auto point
                    = this->layout.fieldNodeCoordinates(Bx, Point<double, 2>{0., 0.}, ix, iy);

                Bx(ix, iy) = std::tanh(0.2 * point[0]) * std::sinh(0.2 * point[1]);
                Jy(ix, iy) = std::sinh(0.3 * point[0]) * std::sinh(0.3 * point[1]);
            }
        }

        for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
        {
            for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
            {
                auto point
                    = this->layout.fieldNodeCoordinates(By, Point<double, 2>{0., 0.}, ix, iy);

                By(ix, iy) = std::tanh(0.4 * point[0]) * std::sinh(0.4 * point[1]);
                Jx(ix, iy) = std::sinh(0.1 * point[0]) * std::sinh(0.1 * point[1]);
            }
        }

        for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
        {
            for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
            {
                auto point
                    = this->layout.fieldNodeCoordinates(Bz, Point<double, 2>{0., 0.}, ix, iy);

                Bz(ix, iy) = std::tanh(0.3 * point[0]) * std::sinh(0.3 * point[1]);
            }
        }
    }
};



class Ohm3DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<3, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> n;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Vx;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Vy;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Vz;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> P;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Jz;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Exnew;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Eynew;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Eznew;
    VecField<NdArrayVector<3>, HybridQuantity> V;
    VecField<NdArrayVector<3>, HybridQuantity> B;
    VecField<NdArrayVector<3>, HybridQuantity> J;
    VecField<NdArrayVector<3>, HybridQuantity> Enew;
    Ohm<GridLayout<GridLayoutImpl>> ohm;

public:
    Ohm3DTest()
        : layout{{{0.1, 0.2, 0.3}}, {{50, 30, 40}}, Point<double, 3>{0., 0., 0.}}
        , n{"n", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho)}
        , Vx{"Vx", HybridQuantity::Scalar::Vx, layout.allocSize(HybridQuantity::Scalar::Vx)}
        , Vy{"Vy", HybridQuantity::Scalar::Vy, layout.allocSize(HybridQuantity::Scalar::Vy)}
        , Vz{"Vz", HybridQuantity::Scalar::Vz, layout.allocSize(HybridQuantity::Scalar::Vz)}
        , P{"P", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , Exnew{"Exnew", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Eynew{"Eynew", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Eznew{"Eznew", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
        , V{"V", HybridQuantity::Vector::V}
        , B{"B", HybridQuantity::Vector::B}
        , J{"J", HybridQuantity::Vector::J}
        , Enew{"Enew", HybridQuantity::Vector::E}
    {
        V.setBuffer("V_x", &Vx);
        V.setBuffer("V_y", &Vy);
        V.setBuffer("V_z", &Vz);
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        J.setBuffer("J_x", &Jx);
        J.setBuffer("J_y", &Jy);
        J.setBuffer("J_z", &Jz);
        Enew.setBuffer("Enew_x", &Exnew);
        Enew.setBuffer("Enew_y", &Eynew);
        Enew.setBuffer("Enew_z", &Eznew);

        auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
        auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
        auto gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
        auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
        auto gsi_p_Z = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Z);
        auto gei_p_Z = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Z);

        auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
        auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
        auto gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
        auto gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);
        auto gsi_d_Z = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Z);
        auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

        for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
        {
            for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
            {
                for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                {
                    auto point = this->layout.fieldNodeCoordinates(
                        Vx, Point<double, 3>{0., 0., 0.}, ix, iy,
                        iz); // only 1 point as all quantities are primal,primal

                    Vx(ix, iy, iz) = std::cosh(0.2 * point[0]) * std::sinh(0.2 * point[1])
                                     * std::tanh(0.2 * point[2]);
                    Vy(ix, iy, iz) = std::cosh(0.4 * point[0]) * std::sinh(0.4 * point[1])
                                     * std::tanh(0.4 * point[2]);
                    Vz(ix, iy, iz) = std::cosh(0.3 * point[0]) * std::sinh(0.3 * point[1])
                                     * std::tanh(0.3 * point[2]);
                    n(ix, iy, iz) = std::exp(-0.1 * point[0]) * std::exp(-0.1 * point[1])
                                    * std::exp(-0.1 * point[2]);
                    P(ix, iy, iz) = std::exp(-0.2 * point[0]) * std::exp(-0.2 * point[1])
                                    * std::exp(-0.2 * point[2]);
                }
                for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
                {
                    auto point = this->layout.fieldNodeCoordinates(Jz, Point<double, 3>{0., 0., 0.},
                                                                   ix, iy, iz);

                    Jz(ix, iy, iz) = std::sinh(0.2 * point[0]) * std::sinh(0.2 * point[1])
                                     * std::tanh(0.2 * point[2]);
                }
            }

            for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
            {
                for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                {
                    auto point = this->layout.fieldNodeCoordinates(Jy, Point<double, 3>{0., 0., 0.},
                                                                   ix, iy, iz);

                    Jy(ix, iy, iz) = std::sinh(0.3 * point[0]) * std::sinh(0.3 * point[1])
                                     * std::tanh(0.3 * point[2]);
                }

                for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
                {
                    auto point = this->layout.fieldNodeCoordinates(Bx, Point<double, 3>{0., 0., 0.},
                                                                   ix, iy, iz);

                    Bx(ix, iy, iz) = std::tanh(0.2 * point[0]) * std::sinh(0.2 * point[1])
                                     * std::cosh(0.2 * point[2]);
                }
            }
        }

        for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
        {
            for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
            {
                for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                {
                    auto point = this->layout.fieldNodeCoordinates(Jx, Point<double, 3>{0., 0., 0.},
                                                                   ix, iy, iz);

                    Jx(ix, iy, iz) = std::sinh(0.1 * point[0]) * std::sinh(0.1 * point[1])
                                     * std::tanh(0.1 * point[2]);
                }
                for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
                {
                    auto point = this->layout.fieldNodeCoordinates(By, Point<double, 3>{0., 0., 0.},
                                                                   ix, iy, iz);

                    By(ix, iy, iz) = std::tanh(0.4 * point[0]) * std::sinh(0.4 * point[1])
                                     * std::cosh(0.4 * point[2]);
                }
            }

            for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
            {
                for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                {
                    auto point = this->layout.fieldNodeCoordinates(Bz, Point<double, 3>{0., 0., 0.},
                                                                   ix, iy, iz);

                    Bz(ix, iy, iz) = std::tanh(0.3 * point[0]) * std::sinh(0.3 * point[1])
                                     * std::cosh(0.3 * point[2]);
                }
            }
        }
    }
};



TEST_F(Ohm1DTest, Ohm1DCalculatedOk)
{
    auto filename_ohmx = std::string{"ohmx_yee_1D_order1.txt"};
    auto filename_ohmy = std::string{"ohmy_yee_1D_order1.txt"};
    auto filename_ohmz = std::string{"ohmz_yee_1D_order1.txt"};
    auto expected_ohmx = read(filename_ohmx);
    auto expected_ohmy = read(filename_ohmy);
    auto expected_ohmz = read(filename_ohmz);

    ohm.setLayout(&layout);
    ohm(n, V, P, B, J, Enew);

    auto psi_X = this->layout.physicalStartIndex(Exnew, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(Exnew, Direction::X);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        EXPECT_THAT(Exnew(ix), ::testing::DoubleNear((expected_ohmx[ix]), 1e-12));
    }

    psi_X = this->layout.physicalStartIndex(Eynew, Direction::X);
    pei_X = this->layout.physicalEndIndex(Eynew, Direction::X);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        EXPECT_THAT(Eynew(ix), ::testing::DoubleNear((expected_ohmy[ix]), 1e-12));
    }

    psi_X = this->layout.physicalStartIndex(Eznew, Direction::X);
    pei_X = this->layout.physicalEndIndex(Eznew, Direction::X);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        EXPECT_THAT(Eznew(ix), ::testing::DoubleNear((expected_ohmz[ix]), 1e-12));
    }
}



TEST_F(Ohm2DTest, Ohm2DCalculatedOk)
{
    auto filename_ohmx = std::string{"ohmx_yee_2D_order1.txt"};
    auto filename_ohmy = std::string{"ohmy_yee_2D_order1.txt"};
    auto filename_ohmz = std::string{"ohmz_yee_2D_order1.txt"};

    auto expected_ohmx = read(filename_ohmx);
    auto expected_ohmy = read(filename_ohmy);
    auto expected_ohmz = read(filename_ohmz);

    ohm.setLayout(&layout);
    ohm(n, V, P, B, J, Enew);

    auto psi_X = this->layout.physicalStartIndex(Exnew, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(Exnew, Direction::X);
    auto psi_Y = this->layout.physicalStartIndex(Exnew, Direction::Y);
    auto pei_Y = this->layout.physicalEndIndex(Exnew, Direction::Y);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ex);
            auto index_ = ix * nPts_[1] + iy;
            // EXPECT_THAT(Exnew(ix, iy), ::testing::DoubleNear((expected_ohmx[index_]), 1e-12));
            EXPECT_FLOAT_EQ(Exnew(ix, iy), expected_ohmx[index_]);
        }
    }

    psi_X = this->layout.physicalStartIndex(Eynew, Direction::X);
    pei_X = this->layout.physicalEndIndex(Eynew, Direction::X);
    psi_Y = this->layout.physicalStartIndex(Eynew, Direction::Y);
    pei_Y = this->layout.physicalEndIndex(Eynew, Direction::Y);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            auto nPts_                   = this->layout.allocSize(HybridQuantity::Scalar::Ey);
            [[maybe_unused]] auto index_ = ix * nPts_[1] + iy;
            // EXPECT_THAT(Eynew(ix, iy), ::testing::DoubleNear((expected_ohmy[index_]), 1e-12));
        }
    }

    psi_X = this->layout.physicalStartIndex(Eznew, Direction::X);
    pei_X = this->layout.physicalEndIndex(Eznew, Direction::X);
    psi_Y = this->layout.physicalStartIndex(Eznew, Direction::Y);
    pei_Y = this->layout.physicalEndIndex(Eznew, Direction::Y);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            auto nPts_                   = this->layout.allocSize(HybridQuantity::Scalar::Ez);
            [[maybe_unused]] auto index_ = ix * nPts_[1] + iy;
            // EXPECT_THAT(Eznew(ix, iy), ::testing::DoubleNear((expected_ohmz[index_]), 1e-12));
        }
    }
}




TEST_F(Ohm3DTest, Ohm3DCalculatedOk)
{
    auto filename_ohmx = std::string{"ohmx_yee_3D_order1.txt"};
    auto filename_ohmy = std::string{"ohmy_yee_3D_order1.txt"};
    auto filename_ohmz = std::string{"ohmz_yee_3D_order1.txt"};

    auto expected_ohmx = read(filename_ohmx);
    auto expected_ohmy = read(filename_ohmy);
    auto expected_ohmz = read(filename_ohmz);

    ohm.setLayout(&layout);
    ohm(n, V, P, B, J, Enew);

    auto psi_X = this->layout.physicalStartIndex(Exnew, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(Exnew, Direction::X);
    auto psi_Y = this->layout.physicalStartIndex(Exnew, Direction::Y);
    auto pei_Y = this->layout.physicalEndIndex(Exnew, Direction::Y);
    auto psi_Z = this->layout.physicalStartIndex(Exnew, Direction::Z);
    auto pei_Z = this->layout.physicalEndIndex(Exnew, Direction::Z);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            for (auto iz = psi_Z; iz <= pei_Z; ++iz)
            {
                auto nPts_                   = this->layout.allocSize(HybridQuantity::Scalar::Ex);
                [[maybe_unused]] auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                // EXPECT_THAT(Exnew(ix, iy, iz),
                //             ::testing::DoubleNear((expected_ohmx[index_]), 1e-12));
            }
        }
    }

    psi_X = this->layout.physicalStartIndex(Eynew, Direction::X);
    pei_X = this->layout.physicalEndIndex(Eynew, Direction::X);
    psi_Y = this->layout.physicalStartIndex(Eynew, Direction::Y);
    pei_Y = this->layout.physicalEndIndex(Eynew, Direction::Y);
    psi_Z = this->layout.physicalStartIndex(Eynew, Direction::Z);
    pei_Z = this->layout.physicalEndIndex(Eynew, Direction::Z);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            for (auto iz = psi_Z; iz <= pei_Z; ++iz)
            {
                auto nPts_                   = this->layout.allocSize(HybridQuantity::Scalar::Ey);
                [[maybe_unused]] auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                // EXPECT_THAT(Eynew(ix, iy, iz),
                //             ::testing::DoubleNear((expected_ohmy[index_]), 1e-12));
            }
        }
    }

    psi_X = this->layout.physicalStartIndex(Eznew, Direction::X);
    pei_X = this->layout.physicalEndIndex(Eznew, Direction::X);
    psi_Y = this->layout.physicalStartIndex(Eznew, Direction::Y);
    pei_Y = this->layout.physicalEndIndex(Eznew, Direction::Y);
    psi_Z = this->layout.physicalStartIndex(Eznew, Direction::Z);
    pei_Z = this->layout.physicalEndIndex(Eznew, Direction::Z);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            for (auto iz = psi_Z; iz <= pei_Z; ++iz)
            {
                auto nPts_                   = this->layout.allocSize(HybridQuantity::Scalar::Ez);
                [[maybe_unused]] auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                // EXPECT_THAT(Eznew(ix, iy, iz),
                //             ::testing::DoubleNear((expected_ohmz[index_]), 1e-12));
            }
        }
    }
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
