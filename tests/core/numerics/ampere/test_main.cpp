#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>
#include <fstream>

#include "phare_core.hpp"

#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/index/index.hpp"

#include "tests/core/data/field/test_field.hpp"
#include "tests/core/data/vecfield/test_vecfield.hpp"
#include "tests/core/data/gridlayout/gridlayout_test.hpp"


using namespace PHARE::core;



struct GridLayoutMock1D
{
    static const auto dimension = 1u;

    template<auto direction>
    double deriv(FieldMock<1> const& /*f*/, MeshIndex<1u> /*mi*/)
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

    template<auto direction>
    double deriv(FieldMock<dimension> const& /*f*/, MeshIndex<2u> /*mi*/)
    {
        return 0;
    }

    std::size_t physicalStartIndex(FieldMock<dimension>&, Direction /*dir*/) { return 0; }

    std::size_t physicalEndIndex(FieldMock<dimension>&, Direction /*dir*/) { return 0; }
};

struct GridLayoutMock3D
{
    static const auto dimension = 3u;

    template<auto direction>
    double deriv(FieldMock<dimension> const& /*f*/, MeshIndex<3u> /*mi*/)
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
    std::size_t constexpr interp = 1;

    auto constexpr check = [](auto ic) {
        auto constexpr dim = ic();

        VecFieldMock<FieldMock<dim>> B, J;
        using GridLayout = GridLayout<GridLayoutImplYee<dim, interp>>;
        Ampere<GridLayout> ampere;
        auto layout = std::make_unique<TestGridLayout<GridLayout>>();
        EXPECT_ANY_THROW(ampere(B, J));
        ampere.setLayout(layout.get());
    };

    check(std::integral_constant<std::size_t, 1>{});
    check(std::integral_constant<std::size_t, 2>{});
    check(std::integral_constant<std::size_t, 3>{});
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
    static constexpr auto dimension    = 1;
    static constexpr auto interp_order = 1;

    using PHARE_TYPES    = PHARE::core::PHARE_Types<dimension, interp_order>;
    using GridLayoutImpl = GridLayoutImplYee<dimension, interp_order>;
    using NdArray_t      = typename PHARE_TYPES::Array_t;

    GridLayout<GridLayoutImpl> layout;

    Field<NdArray_t, HybridQuantity::Scalar> Bx;
    Field<NdArray_t, HybridQuantity::Scalar> By;
    Field<NdArray_t, HybridQuantity::Scalar> Bz;
    Field<NdArray_t, HybridQuantity::Scalar> Jx;
    Field<NdArray_t, HybridQuantity::Scalar> Jy;
    Field<NdArray_t, HybridQuantity::Scalar> Jz;

    VecField<NdArray_t, HybridQuantity> B;
    VecField<NdArray_t, HybridQuantity> J;

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
    static constexpr auto dimension    = 2;
    static constexpr auto interp_order = 1;

    using PHARE_TYPES    = PHARE::core::PHARE_Types<dimension, interp_order>;
    using GridLayoutImpl = GridLayoutImplYee<dimension, interp_order>;
    using NdArray_t      = typename PHARE_TYPES::Array_t;

    GridLayout<GridLayoutImpl> layout;

    Field<NdArray_t, HybridQuantity::Scalar> Bx;
    Field<NdArray_t, HybridQuantity::Scalar> By;
    Field<NdArray_t, HybridQuantity::Scalar> Bz;
    Field<NdArray_t, HybridQuantity::Scalar> Jx;
    Field<NdArray_t, HybridQuantity::Scalar> Jy;
    Field<NdArray_t, HybridQuantity::Scalar> Jz;

    VecField<NdArray_t, HybridQuantity> B;
    VecField<NdArray_t, HybridQuantity> J;

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
    static constexpr auto dimension    = 3;
    static constexpr auto interp_order = 1;

    using PHARE_TYPES    = PHARE::core::PHARE_Types<dimension, interp_order>;
    using GridLayoutImpl = GridLayoutImplYee<dimension, interp_order>;
    using NdArray_t      = typename PHARE_TYPES::Array_t;

    GridLayout<GridLayoutImpl> layout;

    Field<NdArray_t, HybridQuantity::Scalar> Bx;
    Field<NdArray_t, HybridQuantity::Scalar> By;
    Field<NdArray_t, HybridQuantity::Scalar> Bz;
    Field<NdArray_t, HybridQuantity::Scalar> Jx;
    Field<NdArray_t, HybridQuantity::Scalar> Jy;
    Field<NdArray_t, HybridQuantity::Scalar> Jz;

    VecField<NdArray_t, HybridQuantity> B;
    VecField<NdArray_t, HybridQuantity> J;

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
