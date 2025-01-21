#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/point/point.hpp"

#include "tests/core/data/field/test_field.hpp"
#include "tests/core/data/vecfield/test_vecfield.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
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

    std::size_t physicalStartIndex(FieldMock<1>&, Direction /*dir*/) { return 0; }
    std::size_t physicalEndIndex(FieldMock<1>&, Direction /*dir*/) { return 0; }
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

    std::size_t physicalStartIndex(FieldMock<dimension>&, Direction /*dir*/) { return 0; }
    std::size_t physicalEndIndex(FieldMock<dimension>&, Direction /*dir*/) { return 0; }
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
    {
        using GridLayout = GridLayout<GridLayoutImplYee<1, 1>>;
        VecFieldMock<FieldMock<1>> B_1, E_1, Bnew_1;
        Faraday<GridLayout> faraday1d;
        auto layout1d = std::make_unique<TestGridLayout<GridLayout>>();
        EXPECT_ANY_THROW(faraday1d(B_1, E_1, Bnew_1, 1.));
        faraday1d.setLayout(layout1d.get());
    }

    {
        using GridLayout = GridLayout<GridLayoutImplYee<2, 1>>;
        VecFieldMock<FieldMock<2>> B_2, E_2, Bnew_2;
        Faraday<GridLayout> faraday2d;
        auto layout2d = std::make_unique<TestGridLayout<GridLayout>>();
        EXPECT_ANY_THROW(faraday2d(B_2, E_2, Bnew_2, 1.));
        faraday2d.setLayout(layout2d.get());
    }

    {
        using GridLayout = GridLayout<GridLayoutImplYee<3, 1>>;
        VecFieldMock<FieldMock<3>> B_3, E_3, Bnew_3;
        Faraday<GridLayout> faraday3d;
        auto layout3d = std::make_unique<TestGridLayout<GridLayout>>();
        EXPECT_ANY_THROW(faraday3d(B_3, E_3, Bnew_3, 1.));
        faraday3d.setLayout(layout3d.get());
    }
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
    static constexpr auto dim          = 1;
    static constexpr auto interp_order = 1;

    using UsableVecFieldND = UsableVecField<dim>;
    using GridLayoutImpl   = GridLayoutImplYee<dim, interp_order>;

    GridLayout<GridLayoutImpl> layout;

    UsableVecFieldND B, E, Bnew;

    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday1DTest()
        : layout{{{0.1}}, {{50}}, Point{0.}}
        , B{"B", layout, HybridQuantity::Vector::B}
        , E{"E", layout, HybridQuantity::Vector::E}
        , Bnew{"Bnew", layout, HybridQuantity::Vector::B}
    {
    }
};




class Faraday2DTest : public ::testing::Test
{
protected:
    static constexpr auto dim          = 2;
    static constexpr auto interp_order = 1;

    using UsableVecFieldND = UsableVecField<dim>;
    using GridLayoutImpl   = GridLayoutImplYee<dim, interp_order>;

    GridLayout<GridLayoutImpl> layout;

    UsableVecFieldND B, E, Bnew;

    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday2DTest()
        : layout{{{0.1, 0.2}}, {{50, 30}}, Point{0., 0.}}
        , B{"B", layout, HybridQuantity::Vector::B}
        , E{"E", layout, HybridQuantity::Vector::E}
        , Bnew{"Bnew", layout, HybridQuantity::Vector::B}
    {
    }
};




class Faraday3DTest : public ::testing::Test
{
protected:
    static constexpr auto dim          = 3;
    static constexpr auto interp_order = 1;

    using UsableVecFieldND = UsableVecField<dim>;
    using GridLayoutImpl   = GridLayoutImplYee<dim, interp_order>;

    GridLayout<GridLayoutImpl> layout;

    UsableVecFieldND B, E, Bnew;

    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday3DTest()
        : layout{{{0.1, 0.2, 0.3}}, {{50, 30, 40}}, Point{0., 0., 0.}}
        , B{"B", layout, HybridQuantity::Vector::B}
        , E{"E", layout, HybridQuantity::Vector::E}
        , Bnew{"Bnew", layout, HybridQuantity::Vector::B}
    {
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

    auto const& [Ex, Ey, Ez]          = E();
    auto const& [Bx, By, Bz]          = B();
    auto const& [Bxnew, Bynew, Bznew] = Bnew();

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
        if constexpr (std::is_same_v<floater_t<4>, double>)
        {
            EXPECT_THAT(Bynew(ix), ::testing::DoubleNear((expected_dbydt[ix]), 1e-12));
            EXPECT_THAT(Bznew(ix), ::testing::DoubleNear((expected_dbzdt[ix]), 1e-12));
        }
        else
        {
            EXPECT_THAT(Bynew(ix), ::testing::FloatNear((expected_dbydt[ix]), 1e-5));
            EXPECT_THAT(Bznew(ix), ::testing::FloatNear((expected_dbzdt[ix]), 1e-5));
        }
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

    auto const& [Ex, Ey, Ez]          = E();
    auto const& [Bx, By, Bz]          = B();
    auto const& [Bxnew, Bynew, Bznew] = Bnew();

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
            if constexpr (std::is_same_v<floater_t<4>, double>)
            {
                EXPECT_THAT(Bxnew(ix, iy), ::testing::DoubleNear((expected_dbxdt[index_]), 1e-12));
            }
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
            if constexpr (std::is_same_v<floater_t<4>, double>)
            {
                EXPECT_THAT(Bynew(ix, iy), ::testing::DoubleNear((expected_dbydt[index_]), 1e-12));
            }
            else
            {
                // todo
            }
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Bz);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            auto index_ = ix * nPts_[1] + iy;
            if constexpr (std::is_same_v<floater_t<4>, double>)
            {
                EXPECT_THAT(Bznew(ix, iy), ::testing::DoubleNear((expected_dbzdt[index_]), 1e-12));
            }
            else
            {
                // todo
            }
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

    auto const& [Ex, Ey, Ez]          = E();
    auto const& [Bx, By, Bz]          = B();
    auto const& [Bxnew, Bynew, Bznew] = Bnew();

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
                if constexpr (std::is_same_v<floater_t<4>, double>)
                {
                    EXPECT_THAT(Bxnew(ix, iy, iz),
                                ::testing::DoubleNear((expected_dbxdt[index_]), 1e-12));
                }
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
                if constexpr (std::is_same_v<floater_t<4>, double>)
                {
                    EXPECT_THAT(Bynew(ix, iy, iz),
                                ::testing::DoubleNear((expected_dbydt[index_]), 1e-12));
                }
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
                if constexpr (std::is_same_v<floater_t<4>, double>)
                {
                    EXPECT_THAT(Bznew(ix, iy, iz),
                                ::testing::DoubleNear((expected_dbzdt[index_]), 1e-12));
                }
            }
        }
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
