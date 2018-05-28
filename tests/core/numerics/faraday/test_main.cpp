#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "core/numerics/faraday/faraday.h"
#include "data/field/field.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/grid/gridlayoutdefs.h"
#include "data/vecfield/vecfield.h"
#include "utilities/index/index.h"

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
};

struct GridLayoutMock3D
{
    static const int dimension = 3;
};



TEST(Faraday, canBe1D)
{
    Faraday<GridLayoutMock1D> faraday;
}


TEST(Faraday, shouldBeGivenAGridLayoutPointerToBeOperational)
{
    Faraday<GridLayoutMock1D> faraday;
    auto layout = std::make_unique<GridLayoutMock1D>();
    VecFieldMock B, E, Bnew;
    EXPECT_ANY_THROW(faraday(B, E, Bnew));
    faraday.setLayout(layout.get());
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
    uint32 nbrCellsX     = 50;
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
    VecField<NdArrayVector1D<>, HybridQuantity> Bnew;
    VecField<NdArrayVector1D<>, HybridQuantity> E;
    Faraday<GridLayout<GridLayoutImpl>> faraday;

public:
    Faraday1DTest()
        : layout{{{0.1}}, {{nbrCellsX}}, Point<double, 1>{0.}}
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
        , Bnew{"Bnew", HybridQuantity::Vector::B}
        , E{"E", HybridQuantity::Vector::E}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        Bnew.setBuffer("Bnew_x", &Bxnew);
        Bnew.setBuffer("Bnew_y", &Bynew);
        Bnew.setBuffer("Bnew_z", &Bznew);
        E.setBuffer("E_x", &Ex);
        E.setBuffer("E_y", &Ey);
        E.setBuffer("E_z", &Ez);
    }
};




TEST_F(Faraday1DTest, Faraday1DCalculatedOk)
{
    auto filename_dbydt = std::string{"dbydt_yee_1D_order1.txt"};
    auto filename_dbzdt = std::string{"dbzdt_yee_1D_order1.txt"};
    auto expected_dbydt = read(filename_dbydt);
    auto expected_dbzdt = read(filename_dbzdt);

    auto gei_p = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);

    for (auto ix = 0u; ix <= gei_p; ++ix)
    {
        auto point = this->layout.fieldNodeCoordinates(Ey, Point<double, 1>{0.}, ix);
        Ey(ix)     = std::cos(2 * M_PI / 5. * point[0]);
        Ez(ix)     = std::sin(2 * M_PI / 5. * point[0]);
    }

    auto gei_d = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    for (auto ix = 0u; ix <= gei_d; ++ix)
    {
        auto point = this->layout.fieldNodeCoordinates(By, Point<double, 1>{0.}, ix);
        By(ix)     = std::tanh(point[0] - 5. / 2.);
        Bz(ix)     = std::tanh(point[0] - 5. / 2.);
    }

    faraday.setLayout(&layout);
    faraday(B, E, Bnew);

    auto psi_p = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    for (auto ix = psi_p; ix <= pei_p; ++ix)
    {
        EXPECT_THAT(Bynew(ix), ::testing::DoubleNear((expected_dbydt[ix]), 1e-12));
        EXPECT_THAT(Bznew(ix), ::testing::DoubleNear((expected_dbzdt[ix]), 1e-12));
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
