#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "core/data/field/field.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield.h"
#include "core/numerics/ampere/ampere.h"
#include "core/utilities/box/box.h"
#include "core/utilities/index/index.h"

using namespace PHARE::core;

template<std::size_t dim>
struct FieldMock
{
    static auto constexpr dimension = dim;
    double data;
    double& operator()([[maybe_unused]] std::uint32_t i) { return data; }
    double const& operator()([[maybe_unused]] std::uint32_t i) const { return data; }
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



TEST(Ampere, canBe1D)
{
    Ampere<GridLayoutMock1D> ampere;
}

TEST(Ampere, canBe2D)
{
    Ampere<GridLayoutMock2D> ampere;
}

TEST(Ampere, canBe3D)
{
    Ampere<GridLayoutMock3D> ampere;
}

TEST(Ampere, shouldBeGivenAGridLayoutPointerToBeOperational)
{
    VecFieldMock<FieldMock<1>> B, J;

    Ampere<GridLayoutMock1D> ampere1d;
    auto layout1d = std::make_unique<GridLayoutMock1D>();
    EXPECT_ANY_THROW(ampere1d(B, J));
    ampere1d.setLayout(layout1d.get());

    Ampere<GridLayoutMock2D> ampere2d;
    auto layout2d = std::make_unique<GridLayoutMock2D>();
    // EXPECT_ANY_THROW(ampere2d(B2, J2));
    ampere2d.setLayout(layout2d.get());

    Ampere<GridLayoutMock3D> ampere3d;
    auto layout3d = std::make_unique<GridLayoutMock3D>();
    // EXPECT_ANY_THROW(ampere3d(B3, J3));
    ampere3d.setLayout(layout3d.get());
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



class Ampere1DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<1, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr auto interp_order = GridLayoutImpl::interp_order;

    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jz;

    VecField<NdArrayVector<1>, HybridQuantity> B;
    VecField<NdArrayVector<1>, HybridQuantity> J;

    Ampere<GridLayout<GridLayoutImpl>> ampere;

public:
    Ampere1DTest()
        : layout{{{0.1}}, {{50}}, Point{0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , B{"B", HybridQuantity::Vector::B}
        , J{"J", HybridQuantity::Vector::J}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        J.setBuffer("J_x", &Jx);
        J.setBuffer("J_y", &Jy);
        J.setBuffer("J_z", &Jz);
    }
};


class Ampere2DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<2, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr auto interp_order = GridLayoutImpl::interp_order;

    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Jz;

    VecField<NdArrayVector<2>, HybridQuantity> B;
    VecField<NdArrayVector<2>, HybridQuantity> J;

    Ampere<GridLayout<GridLayoutImpl>> ampere;

public:
    Ampere2DTest()
        : layout{{{0.1, 0.2}}, {{50, 30}}, Point{0., 0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , B{"B", HybridQuantity::Vector::B}
        , J{"J", HybridQuantity::Vector::J}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        J.setBuffer("J_x", &Jx);
        J.setBuffer("J_y", &Jy);
        J.setBuffer("J_z", &Jz);
    }
};


class Ampere3DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<3, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr auto interp_order = GridLayoutImpl::interp_order;

    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Jz;

    VecField<NdArrayVector<3>, HybridQuantity> B;
    VecField<NdArrayVector<3>, HybridQuantity> J;

    Ampere<GridLayout<GridLayoutImpl>> ampere;

public:
    Ampere3DTest()
        : layout{{{0.1, 0.2, 0.3}}, {{50, 30, 40}}, Point{0., 0., 0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , B{"B", HybridQuantity::Vector::B}
        , J{"J", HybridQuantity::Vector::J}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        J.setBuffer("J_x", &Jx);
        J.setBuffer("J_y", &Jy);
        J.setBuffer("J_z", &Jz);
    }
};



TEST_F(Ampere1DTest, ampere1DCalculatedOk)
{
    auto filename_jy = std::string{"jy_yee_1D_order1.txt"};
    auto filename_jz = std::string{"jz_yee_1D_order1.txt"};
    auto expectedJy  = read(filename_jy);
    auto expectedJz  = read(filename_jz);

    std::uint32_t gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
    std::uint32_t gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    for (std::uint32_t ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        auto point = this->layout.fieldNodeCoordinates(By, Point{0.}, ix);
        By(ix)     = std::cos(2 * M_PI / 5. * point[0]);
        Bz(ix)     = std::sin(2 * M_PI / 5. * point[0]);
    }

    ampere.setLayout(&layout);
    ampere(B, J);

    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    for (std::uint32_t ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        EXPECT_THAT(Jy(ix), ::testing::DoubleNear((expectedJy[ix]), 1e-12));
        EXPECT_THAT(Jz(ix), ::testing::DoubleNear((expectedJz[ix]), 1e-12));
    }
}


TEST_F(Ampere2DTest, ampere2DCalculatedOk)
{
    auto filename_jx = std::string{"jx_yee_2D_order1.txt"};
    auto filename_jy = std::string{"jy_yee_2D_order1.txt"};
    auto filename_jz = std::string{"jz_yee_2D_order1.txt"};
    auto expectedJx  = read(filename_jx);
    auto expectedJy  = read(filename_jy);
    auto expectedJz  = read(filename_jz);

    std::uint32_t gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
    std::uint32_t gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    std::uint32_t gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
    std::uint32_t gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    std::uint32_t gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
    std::uint32_t gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    std::uint32_t gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
    std::uint32_t gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);

    for (std::uint32_t ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (std::uint32_t iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Bx, Point{0., 0.}, ix, iy);
            Bx(ix, iy) = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    for (std::uint32_t ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (std::uint32_t iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(By, Point{0., 0.}, ix, iy);
            By(ix, iy) = std::cos(2 * M_PI / 5. * point[0]) * std::tanh(2 * M_PI / 6. * point[1]);
        }
    }

    for (std::uint32_t ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (std::uint32_t iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Bz, Point{0., 0.}, ix, iy);
            Bz(ix, iy) = std::sin(2 * M_PI / 5. * point[0]) * std::tanh(2 * M_PI / 6. * point[1]);
        }
    }

    ampere.setLayout(&layout);
    ampere(B, J);

    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

    std::array<std::uint32_t, 2> nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jx);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            std::uint32_t index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(Jx(ix, iy), ::testing::DoubleNear((expectedJx[index_]), 1e-12));
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jy);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            std::uint32_t index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(Jy(ix, iy), ::testing::DoubleNear((expectedJy[index_]), 1e-12));
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jz);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            std::uint32_t index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(Jz(ix, iy), ::testing::DoubleNear((expectedJz[index_]), 1e-12));
        }
    }
}



TEST_F(Ampere3DTest, ampere3DCalculatedOk)
{
    auto filename_jx = std::string{"jx_yee_3D_order1.txt"};
    auto filename_jy = std::string{"jy_yee_3D_order1.txt"};
    auto filename_jz = std::string{"jz_yee_3D_order1.txt"};
    auto expectedJx  = read(filename_jx);
    auto expectedJy  = read(filename_jy);
    auto expectedJz  = read(filename_jz);

    std::uint32_t gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
    std::uint32_t gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    std::uint32_t gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
    std::uint32_t gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    std::uint32_t gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
    std::uint32_t gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    std::uint32_t gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
    std::uint32_t gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);
    std::uint32_t gsi_p_Z = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Z);
    std::uint32_t gei_p_Z = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Z);
    std::uint32_t gsi_d_Z = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Z);
    std::uint32_t gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    for (std::uint32_t ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (std::uint32_t iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            for (std::uint32_t iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
            {
                Point<double, 3> point = this->layout.fieldNodeCoordinates(
                    Bx, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                Bx(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])
                                 * std::cos(2 * M_PI / 6. * point[1])
                                 * std::tanh(2 * M_PI / 12. * point[2]);
            }
        }
    }

    for (std::uint32_t ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (std::uint32_t iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            for (std::uint32_t iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
            {
                Point<double, 3> point = this->layout.fieldNodeCoordinates(
                    By, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                By(ix, iy, iz) = std::tanh(2 * M_PI / 5. * point[0])
                                 * std::sin(2 * M_PI / 6. * point[1])
                                 * std::cos(2 * M_PI / 12. * point[2]);
            }
        }
    }

    for (std::uint32_t ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (std::uint32_t iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            for (std::uint32_t iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
            {
                Point<double, 3> point = this->layout.fieldNodeCoordinates(
                    Bz, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                Bz(ix, iy, iz) = std::cos(2 * M_PI / 5. * point[0])
                                 * std::tanh(2 * M_PI / 6. * point[1])
                                 * std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    ampere.setLayout(&layout);
    ampere(B, J);

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

    std::array<std::uint32_t, 3> nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jx);

    for (std::uint32_t ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (std::uint32_t iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (std::uint32_t iz = psi_p_Z; iz <= pei_p_Z; ++iz)
            {
                std::uint32_t index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(Jx(ix, iy, iz), ::testing::DoubleNear((expectedJx[index_]), 1e-12));
            }
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jy);

    for (std::uint32_t ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (std::uint32_t iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            for (std::uint32_t iz = psi_p_Z; iz <= pei_p_Z; ++iz)
            {
                std::uint32_t index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(Jy(ix, iy, iz), ::testing::DoubleNear((expectedJy[index_]), 1e-12));
            }
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jz);

    for (std::uint32_t ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (std::uint32_t iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (std::uint32_t iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                std::uint32_t index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(Jz(ix, iy, iz), ::testing::DoubleNear((expectedJz[index_]), 1e-12));
            }
        }
    }
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
