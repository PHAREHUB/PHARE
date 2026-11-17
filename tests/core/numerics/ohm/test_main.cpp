
#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/pressure_term/pressure_term.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"


#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <variant>

using namespace PHARE::core;


template<int dim, int interp>
class NDlayout
{
    NDlayout() {}

    using nDL = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;

public:
    static nDL create()
    {
        if constexpr (dim == 1)
        {
            return {{{0.1}}, {{50}}, {0.}};
        }
        else if constexpr (dim == 2)
        {
            return {{{0.1, 0.2}}, {{50, 40}}, {0., 0.}};
        }
        else if constexpr (dim == 3)
        {
            return {{{0.1, 0.2, 0.3}}, {{50, 40, 30}}, {0., 0., 0.}};
        }
    }
};




PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;

    dict["resistivity"]       = 1.0;
    dict["hyper_resistivity"] = 0.01;

    return dict;
}



template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct OhmTest : public ::testing::Test
{
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();

    using GridYee = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;
    using UsableVecFieldND = UsableVecField<dim>;
    using Grid_t           = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using Ohm_t            = Ohm<GridYee>;

    GridYee layout = NDlayout<dim, interp>::create();

    // electron temperature of the isothermal pressure term, as in test_ohm.py
    static constexpr double Te = 0.09;

    Grid_t n;
    UsableVecFieldND V, B, J, Enew, pressureTerm;
    UsableVecFieldND Ji; // ion current on the E nodes
    UsableVecFieldND Ge; // electron flux Ji - J on the E nodes, the input of Ohm's ideal term

    // what Electrons::computeFlux does, on the whole ghost box: G = Ji - J
    auto& electronFlux()
    {
        for (auto const c : {Component::X, Component::Y, Component::Z})
            layout.evalOnGhostBox(
                Ge(c), [&](auto const&... ijk) { Ge(c)(ijk...) = Ji(c)(ijk...) - J(c)(ijk...); });
        return Ge;
    }

    OhmInfo ohm_info;

    OhmTest()
        : n{"n", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho)}
        , V{"V", layout, HybridQuantity::Vector::V}
        , B{"B", layout, HybridQuantity::Vector::B}
        , J{"J", layout, HybridQuantity::Vector::J}
        , Enew{"Enew", layout, HybridQuantity::Vector::E}
        , pressureTerm{"pressureTerm", layout, HybridQuantity::Vector::E}
        , Ji{"Ji", layout, HybridQuantity::Vector::J}
        , Ge{"Ge", layout, HybridQuantity::Vector::J}
        , ohm_info{OhmInfo::FROM(createDict())}
    {
        auto const& [Bx, By, Bz]          = B();
        auto const& [Jx, Jy, Jz]          = J();
        auto const& [Vx, Vy, Vz]          = V();
        auto const& [Exnew, Eynew, Eznew] = Enew();

        if constexpr (dim == 1)
        {
            auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
            auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
            auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
            auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

            for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
            {
                auto point = this->layout.fieldNodeCoordinates(
                    n, this->layout.localToAMR(Point{ix}.as_signed()));

                n(ix)  = std::cosh(0.5 * point[0]);
                Vx(ix) = std::sinh(0.2 * point[0]);
                Vy(ix) = std::sinh(0.3 * point[0]);
                Vz(ix) = std::sinh(0.4 * point[0]);
                Bx(ix) = std::cosh(0.2 * point[0]);
                Jy(ix) = std::tanh(0.3 * point[0]);
                Jz(ix) = std::tanh(0.4 * point[0]);
            }

            for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
            {
                auto point = this->layout.fieldNodeCoordinates(
                    Bz, this->layout.localToAMR(Point{ix}.as_signed()));

                By(ix) = std::cosh(0.3 * point[0]);
                Bz(ix) = std::cosh(0.4 * point[0]);
                Jx(ix) = std::tanh(0.2 * point[0]);
            }
        }

        if constexpr (dim == 2)
        {
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
                        n, this->layout.localToAMR(Point{ix, iy}.as_signed()));

                    n(ix, iy)  = std::cosh(0.5 * point[0]) * std::cosh(0.5 * point[1]);
                    Vx(ix, iy) = std::sinh(0.2 * point[0]) * std::sinh(0.2 * point[1]);
                    Vy(ix, iy) = std::sinh(0.3 * point[0]) * std::sinh(0.3 * point[1]);
                    Vz(ix, iy) = std::sinh(0.4 * point[0]) * std::sinh(0.4 * point[1]);
                    Jz(ix, iy) = std::tanh(0.4 * point[0]) * std::tanh(0.4 * point[1]);
                }
                for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
                {
                    auto point = this->layout.fieldNodeCoordinates(
                        Bx, this->layout.localToAMR(Point{ix, iy}.as_signed()));

                    Bx(ix, iy) = std::cosh(0.2 * point[0]) * std::cosh(0.2 * point[1]);
                    Jy(ix, iy) = std::tanh(0.3 * point[0]) * std::tanh(0.3 * point[1]);
                }
            }
            for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
            {
                for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
                {
                    auto point = this->layout.fieldNodeCoordinates(
                        Jx, this->layout.localToAMR(Point{ix, iy}.as_signed()));

                    By(ix, iy) = std::cosh(0.3 * point[0]) * std::cosh(0.3 * point[1]);
                    Jx(ix, iy) = std::tanh(0.2 * point[0]) * std::tanh(0.2 * point[1]);
                }
                for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
                {
                    auto point = this->layout.fieldNodeCoordinates(
                        Bz, this->layout.localToAMR(Point{ix, iy}.as_signed()));

                    Bz(ix, iy) = std::cosh(0.4 * point[0]) * std::cosh(0.4 * point[1]);
                }
            }
        }

        if constexpr (dim == 3)
        {
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
                            n, this->layout.localToAMR(Point{ix, iy, iz}.as_signed()));

                        n(ix, iy, iz) = std::cosh(0.5 * point[0]) * std::cosh(0.5 * point[1])
                                        * std::cosh(0.5 * point[2]);
                        Vx(ix, iy, iz) = std::sinh(0.2 * point[0]) * std::sinh(0.2 * point[1])
                                         * std::sinh(0.2 * point[2]);
                        Vy(ix, iy, iz) = std::sinh(0.3 * point[0]) * std::sinh(0.3 * point[1])
                                         * std::sinh(0.3 * point[2]);
                        Vz(ix, iy, iz) = std::sinh(0.4 * point[0]) * std::sinh(0.4 * point[1])
                                         * std::sinh(0.4 * point[2]);
                    }
                    for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            Jz, this->layout.localToAMR(Point{ix, iy, iz}.as_signed()));

                        Jz(ix, iy, iz) = std::tanh(0.4 * point[0]) * std::tanh(0.4 * point[1])
                                         * std::tanh(0.4 * point[2]);
                    }
                }
                for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
                {
                    for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            Jy, this->layout.localToAMR(Point{ix, iy, iz}.as_signed()));

                        Jy(ix, iy, iz) = std::tanh(0.3 * point[0]) * std::tanh(0.3 * point[1])
                                         * std::tanh(0.3 * point[2]);
                    }
                    for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            Bx, this->layout.localToAMR(Point{ix, iy, iz}.as_signed()));

                        Bx(ix, iy, iz) = std::cosh(0.2 * point[0]) * std::cosh(0.2 * point[1])
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
                        auto point = this->layout.fieldNodeCoordinates(
                            Jx, this->layout.localToAMR(Point{ix, iy, iz}.as_signed()));

                        Jx(ix, iy, iz) = std::tanh(0.2 * point[0]) * std::tanh(0.2 * point[1])
                                         * std::tanh(0.2 * point[2]);
                    }
                    for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            By, this->layout.localToAMR(Point{ix, iy, iz}.as_signed()));

                        By(ix, iy, iz) = std::cosh(0.3 * point[0]) * std::cosh(0.3 * point[1])
                                         * std::cosh(0.3 * point[2]);
                    }
                }
                for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
                {
                    for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            Bz, this->layout.localToAMR(Point{ix, iy, iz}.as_signed()));

                        Bz(ix, iy, iz) = std::cosh(0.4 * point[0]) * std::cosh(0.4 * point[1])
                                         * std::cosh(0.4 * point[2]);
                    }
                }
            }
        }

        // Ji: the profiles of V, on the E nodes (as test_ohm.py)
        double const k[3] = {0.2, 0.3, 0.4};
        for (std::uint8_t c = 0; c < 3; ++c)
        {
            auto& field = Ji[c];
            layout.evalOnGhostBox(field, [&](auto const&... ijk) {
                auto const point = layout.fieldNodeCoordinates(
                    field, layout.localToAMR(Point{ijk...}.as_signed()));
                double value = 1.;
                for (std::size_t d = 0; d < dim; ++d)
                    value *= std::sinh(k[c] * point[d]);
                field(ijk...) = value;
            });
        }
    }

    ~OhmTest() {}
};



using OhmTupleInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<3>, InterpConst<1>>>;

TYPED_TEST_SUITE(OhmTest, OhmTupleInfos);




TYPED_TEST(OhmTest, ThatOhmInfoCanBeBuiltFromDict)
{
    auto info = OhmInfo::FROM(createDict());
}


// asEnumConstant enumerates [0, Enum::count), which is what replaced the EnumTraits value table
static_assert(std::is_same_v<EnumConstantVariant_t<HyperMode>,
                             std::variant<std::integral_constant<HyperMode, HyperMode::constant>,
                                          std::integral_constant<HyperMode, HyperMode::spatial>>>);

TYPED_TEST(OhmTest, ThatAsEnumConstantLiftsEveryHyperMode)
{
    for (auto const mode : {HyperMode::constant, HyperMode::spatial})
        std::visit([&](auto const tag) { EXPECT_EQ(mode, decltype(tag)::value); },
                   asEnumConstant(mode));

    EXPECT_THROW(asEnumConstant(HyperMode::count), std::runtime_error);
}


// hyper_mode is written to the dict as the int of the HyperMode member exposed to python by
// cpp_etc (see pharein.initialize.add_enum_int), and read back as the enum
TYPED_TEST(OhmTest, ThatHyperModeIsReadFromDictAsAnEnum)
{
    EXPECT_EQ(HyperMode::constant, OhmInfo::FROM(createDict()).hyper_mode); // absent, defaulted

    auto dict          = createDict();
    dict["hyper_mode"] = static_cast<int>(HyperMode::spatial);
    EXPECT_EQ(HyperMode::spatial, OhmInfo::FROM(dict).hyper_mode);

    dict["hyper_mode"] = static_cast<int>(HyperMode::constant);
    EXPECT_EQ(HyperMode::constant, OhmInfo::FROM(dict).hyper_mode);
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



TYPED_TEST(OhmTest, ThatElectricFieldIsOkFromOhmsLaw)
{
    using Ohm_t = TestFixture::Ohm_t;

    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    std::string prefix{"ohm"};

    auto filenameX
        = prefix + "x_yee_" + std::to_string(dim) + "D_order" + std::to_string(interp) + ".txt";
    auto filenameY
        = prefix + "y_yee_" + std::to_string(dim) + "D_order" + std::to_string(interp) + ".txt";
    auto filenameZ
        = prefix + "z_yee_" + std::to_string(dim) + "D_order" + std::to_string(interp) + ".txt";
    auto expected_ohmX = read(filenameX);
    auto expected_ohmY = read(filenameY);
    auto expected_ohmZ = read(filenameZ);

    auto const layout = NDlayout<dim, interp>::create();

    // the reference values (test_ohm.py) use the isothermal pressure term, so build it the same way
    IsothermalPressureTerm<typename TestFixture::GridYee>{TestFixture::Te}(layout, this->n,
                                                                           this->pressureTerm);
    Ohm_t{this->ohm_info, layout}(this->n, this->electronFlux(), this->pressureTerm, this->B,
                                  this->J, this->Enew);


    auto const& [Exnew, Eynew, Eznew] = this->Enew();

    if constexpr (dim == 1)
    {
        auto psi_X = this->layout.physicalStartIndex(Exnew, Direction::X);
        auto pei_X = this->layout.physicalEndIndex(Exnew, Direction::X);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            EXPECT_THAT(Exnew(ix), ::testing::DoubleNear((expected_ohmX[ix]), 1e-12));
        }

        psi_X = this->layout.physicalStartIndex(Eynew, Direction::X);
        pei_X = this->layout.physicalEndIndex(Eynew, Direction::X);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            EXPECT_THAT(Eynew(ix), ::testing::DoubleNear((expected_ohmY[ix]), 1e-12));
        }

        psi_X = this->layout.physicalStartIndex(Eznew, Direction::X);
        pei_X = this->layout.physicalEndIndex(Eznew, Direction::X);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            EXPECT_THAT(Eznew(ix), ::testing::DoubleNear((expected_ohmZ[ix]), 1e-12));
        }
    }

    if constexpr (dim == 2)
    {
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
                EXPECT_THAT(Exnew(ix, iy), ::testing::DoubleNear((expected_ohmX[index_]), 1e-12));
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
                auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ey);
                auto index_ = ix * nPts_[1] + iy;
                EXPECT_THAT(Eynew(ix, iy), ::testing::DoubleNear((expected_ohmY[index_]), 1e-12));
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
                auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ez);
                auto index_ = ix * nPts_[1] + iy;
                EXPECT_THAT(Eznew(ix, iy), ::testing::DoubleNear((expected_ohmZ[index_]), 1e-12));
            }
        }
    }

    if constexpr (dim == 3)
    {
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
                    auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ex);
                    auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                    EXPECT_THAT(Exnew(ix, iy, iz),
                                ::testing::DoubleNear((expected_ohmX[index_]), 1e-10));
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
                    auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ey);
                    auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                    EXPECT_THAT(Eynew(ix, iy, iz),
                                ::testing::DoubleNear((expected_ohmY[index_]), 1e-10));
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
                    auto nPts_  = this->layout.allocSize(HybridQuantity::Scalar::Ez);
                    auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                    EXPECT_THAT(Eznew(ix, iy, iz),
                                ::testing::DoubleNear((expected_ohmZ[index_]), 1e-10));
                }
            }
        }
    }
}


// Ohm adds the pressure term it is given, whatever the closure that computed it: with
// Ji = B = J = 0 and no resistivity, E is the pressure term, value for value
TYPED_TEST(OhmTest, ThatElectricFieldIsThePressureTermWhenNothingElseActs)
{
    using Ohm_t = TestFixture::Ohm_t;

    auto& layout    = this->layout;
    auto const zero = [&](auto& field) {
        layout.evalOnGhostBox(field, [&](auto const&... ijk) { field(ijk...) = 0.; });
    };
    for (auto* vecfield : {&this->Ji, &this->B, &this->J})
        for (auto const component : {Component::X, Component::Y, Component::Z})
            zero((*vecfield)(component));

    auto k = 0;
    for (auto const c : {Component::X, Component::Y, Component::Z})
    {
        auto& term = this->pressureTerm(c);
        layout.evalOnBox(term, [&](auto const&... ijk) {
            auto const hash = std::sin(12.9898 * ++k) * 43758.5453;
            term(ijk...)    = hash - std::floor(hash) - 0.5;
        });
    }

    auto const dict = [] {
        PHARE::initializer::PHAREDict dict;
        dict["resistivity"]       = 0.0;
        dict["hyper_resistivity"] = 0.0;
        return dict;
    }();
    Ohm_t{OhmInfo::FROM(dict), layout}(this->n, this->electronFlux(), this->pressureTerm, this->B,
                                       this->J, this->Enew);

    for (auto const c : {Component::X, Component::Y, Component::Z})
    {
        auto& Ec   = this->Enew(c);
        auto& term = this->pressureTerm(c);
        layout.evalOnBox(Ec, [&](auto const&... ijk) { EXPECT_EQ(Ec(ijk...), term(ijk...)); });
    }
}


// ---------------------------------------------------------------------------------------------
// The ideal term does no work on the electrons (#1577). With G = Ji - J the electron flux on the E
// nodes, E = -G x B / n must give sum_i G_i E_i = 0 over a periodic grid, for any n > 0, Ji, J,
// B. The fields are periodic noise of the AMR index (period = number of cells), so the ghost
// nodes are periodic copies. The former stencil ("Ve = (F - Proj(J))/n on the moment nodes,
// Ve and B averaged to the E node, then multiplied") is run in the same harness, with G = Avg(F)
// - J built by the test's own averaging, and must fail, so the test can tell the two apart.
// ---------------------------------------------------------------------------------------------

template<typename TypeInfo>
struct OhmWorkTest : public ::testing::Test
{
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();

    using GridYee = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;
    using Grid_t  = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using Ohm_t   = Ohm<GridYee>;
    using Scalar  = HybridQuantity::Scalar;

    GridYee layout = NDlayout<dim, interp>::create();
    std::array<int, dim> const period
        = for_N_make_array<dim>([&](auto i) { return int(layout.nbrCells()[i]); });

    Grid_t n{"n", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho)};
    Grid_t Vex{"Vex", Scalar::Vx, layout.allocSize(Scalar::Vx)};
    Grid_t Vey{"Vey", Scalar::Vy, layout.allocSize(Scalar::Vy)};
    Grid_t Vez{"Vez", Scalar::Vz, layout.allocSize(Scalar::Vz)};
    UsableVecField<dim> F{"F", layout, HybridQuantity::Vector::V};
    UsableVecField<dim> Ji{"Ji", layout, HybridQuantity::Vector::J};
    UsableVecField<dim> Ge{"Ge", layout, HybridQuantity::Vector::J}; // Ji - J, Ohm's input
    UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
    UsableVecField<dim> J{"J", layout, HybridQuantity::Vector::J};
    UsableVecField<dim> E{"E", layout, HybridQuantity::Vector::E};
    UsableVecField<dim> pressureTerm{"pressureTerm", layout, HybridQuantity::Vector::E};

    // uniform in [-1, 1), periodic in the AMR index (splitmix64 of the wrapped index)
    double noise(std::uint64_t seed, auto const& amr) const
    {
        std::uint64_t h = seed * 0x9E3779B97F4A7C15ull;
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const wrapped = ((amr[d] % period[d]) + period[d]) % period[d];
            h ^= std::uint64_t(wrapped) + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
            h = (h ^ (h >> 30)) * 0xBF58476D1CE4E5B9ull;
            h = (h ^ (h >> 27)) * 0x94D049BB133111EBull;
            h ^= h >> 31;
        }
        return 2. * double(h >> 11) / double(1ull << 53) - 1.;
    }

    void fill(auto& field, std::uint64_t seed, auto&& fn)
    {
        layout.ghost_amr_lcl_idx(
            field, [&](auto const& amr, auto const& lcl) { field(lcl) = fn(noise(seed, amr)); });
    }

    OhmWorkTest()
    {
        // n log-normal over 0.05 to 20, B and F with nonzero means
        fill(n, 1, [](double u) { return std::exp(3. * u); });
        double const Bmean[3] = {0.3, -0.7, 1.1}, Fmean[3] = {0.2, 0.5, -0.4};
        for (std::uint8_t c = 0; c < 3; ++c)
        {
            fill(B[c], 10 + c, [&](double u) { return Bmean[c] + u; });
            fill(F[c], 20 + c, [&](double u) { return Fmean[c] + u; });
            fill(Ji[c], 60 + c, [&](double u) { return Fmean[c] + u; });
            fill(J[c], 30 + c, [](double u) { return u; });
            fill(pressureTerm[c], 40 + c, [](double) { return 0.; });
            fill(E[c], 50 + c, [](double) { return 0.; });
        }
    }

    // G_c on the E_c nodes: Ji - J, or for the former stencil Avg(F) - J with the test's own
    // averaging: E_c is dual along c (if c < dim) and primal elsewhere, as are the moments;
    // dual node i sits between primal nodes i and i + 1
    bool former = false;
    double G(std::uint8_t c, auto const& lcl) const
    {
        if (!former)
            return Ji[c](lcl) - J[c](lcl);
        auto f = F[c](lcl);
        if (c < dim)
        {
            auto up = lcl;
            ++up[c];
            f = 0.5 * (f + F[c](up));
        }
        return f - J[c](lcl);
    }

    // sum_i G_i E_i over each E component's own periodic nodes (upper primal node excluded),
    // and the sum of |G_i E_i| as the scale
    std::pair<double, double> work() const
    {
        double sum = 0., scale = 0.;
        for (std::uint8_t c = 0; c < 3; ++c)
            layout.domain_amr_lcl_idx(E[c], [&](auto const& amr, auto const& lcl) {
                for (std::size_t d = 0; d < dim; ++d)
                    if (amr[d] == period[d])
                        return;
                auto const ge = G(c, lcl) * E[c](lcl);
                sum += ge;
                scale += std::abs(ge);
            });
        return {sum, scale};
    }

    // the former ideal term: Ve = (F - Proj(J)) / n on the moment nodes, then Ve and B averaged
    // separately to each E node and multiplied there
    void formerIdealTerm()
    {
        former            = true;
        auto const shrink = for_N_make_array<dim>([](auto) { return std::uint32_t{1}; });
        layout.evalOnShrinkedGhostBox(n, shrink, [&](auto const&... args) {
            auto const i = MeshIndex<dim>{args...};
            Vex(i) = (F[0](i) - GridYee::template project<Scalar::Jx, Scalar::Vx>(J[0], i)) / n(i);
            Vey(i) = (F[1](i) - GridYee::template project<Scalar::Jy, Scalar::Vy>(J[1], i)) / n(i);
            Vez(i) = (F[2](i) - GridYee::template project<Scalar::Jz, Scalar::Vz>(J[2], i)) / n(i);
        });
        auto P = [&]<auto src, auto dst>(auto const& f, auto const& i) {
            return GridYee::template project<src, dst>(f, i);
        };
        layout.evalOnBox(E[0], [&](auto const&... args) {
            auto const i = MeshIndex<dim>{args...};
            E[0](i)      = -P.template operator()<Scalar::Vy, Scalar::Ex>(Vey, i)
                          * P.template operator()<Scalar::Bz, Scalar::Ex>(B[2], i)
                      + P.template operator()<Scalar::Vz, Scalar::Ex>(Vez, i)
                            * P.template operator()<Scalar::By, Scalar::Ex>(B[1], i);
        });
        layout.evalOnBox(E[1], [&](auto const&... args) {
            auto const i = MeshIndex<dim>{args...};
            E[1](i)      = -P.template operator()<Scalar::Vz, Scalar::Ey>(Vez, i)
                          * P.template operator()<Scalar::Bx, Scalar::Ey>(B[0], i)
                      + P.template operator()<Scalar::Vx, Scalar::Ey>(Vex, i)
                            * P.template operator()<Scalar::Bz, Scalar::Ey>(B[2], i);
        });
        layout.evalOnBox(E[2], [&](auto const&... args) {
            auto const i = MeshIndex<dim>{args...};
            E[2](i)      = -P.template operator()<Scalar::Vx, Scalar::Ez>(Vex, i)
                          * P.template operator()<Scalar::By, Scalar::Ez>(B[1], i)
                      + P.template operator()<Scalar::Vy, Scalar::Ez>(Vey, i)
                            * P.template operator()<Scalar::Bx, Scalar::Ez>(B[0], i);
        });
    }
};

using OhmWorkTupleInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<1>, InterpConst<2>>,
                     std::pair<DimConst<1>, InterpConst<3>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<2>, InterpConst<2>>, std::pair<DimConst<2>, InterpConst<3>>,
                     std::pair<DimConst<3>, InterpConst<1>>, std::pair<DimConst<3>, InterpConst<2>>,
                     std::pair<DimConst<3>, InterpConst<3>>>;

TYPED_TEST_SUITE(OhmWorkTest, OhmWorkTupleInfos);


TYPED_TEST(OhmWorkTest, ThatTheIdealTermDoesNoWorkOnTheElectrons)
{
    using Ohm_t     = TestFixture::Ohm_t;
    auto const dict = [] {
        PHARE::initializer::PHAREDict dict;
        dict["resistivity"]       = 0.0;
        dict["hyper_resistivity"] = 0.0;
        return dict;
    }();

    for (auto const c : {Component::X, Component::Y, Component::Z})
        this->layout.evalOnGhostBox(this->Ge(c), [&](auto const&... ijk) {
            this->Ge(c)(ijk...) = this->Ji(c)(ijk...) - this->J(c)(ijk...);
        });
    Ohm_t{OhmInfo::FROM(dict), this->layout}(this->n, this->Ge, this->pressureTerm, this->B,
                                             this->J, this->E);
    auto const [sum, scale] = this->work();
    EXPECT_GT(scale, 0.);
    EXPECT_LT(std::abs(sum), 1e-13 * scale) << "relative work " << sum / scale;
    std::cout << "ideal term, relative work " << sum / scale << "\n";
}


TYPED_TEST(OhmWorkTest, ThatTheFormerIdealTermDoesWork)
{
    this->formerIdealTerm();
    auto const [sum, scale] = this->work();
    EXPECT_GT(std::abs(sum), 1e-6 * scale) << "relative work " << sum / scale;
    std::cout << "former ideal term, relative work " << sum / scale << "\n";
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
