#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "data/field/field.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/grid/gridlayoutdefs.h"
#include "data/vecfield/vecfield.h"
#include "numerics/ampere/ampere.h"
#include "utilities/index/index.h"

using namespace PHARE;

struct FieldMock
{
    double data;
    double& operator()(uint32 i) { return data; }
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



TEST(Ampere, canBe1D)
{
    Ampere<GridLayoutMock1D> ampere;
}


TEST(Ampere, shouldBeGivenAGridLayoutPointerToBeOperational)
{
    Ampere<GridLayoutMock1D> ampere;
    auto layout = std::make_unique<GridLayoutMock1D>();
    VecFieldMock B, J;
    EXPECT_ANY_THROW(ampere(B, J));
    ampere.setLayout(layout.get());
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
    uint32 nbrCellsX     = 50;
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Bx;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> By;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Bz;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Jx;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Jy;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Jz;
    VecField<NdArrayVector1D<>, HybridQuantity> B;
    VecField<NdArrayVector1D<>, HybridQuantity> J;
    Ampere<GridLayout<GridLayoutImpl>> ampere;

public:
    Ampere1DTest()
        : layout{{{0.1}}, {{nbrCellsX}}, Point<double, 1>{0.}}
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

    auto gei_d = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    for (auto ix = 0u; ix <= gei_d; ++ix)
    {
        auto point = this->layout.fieldNodeCoordinates(By, Point<double, 1>{0.}, ix);
        By(ix)     = std::cos(2 * M_PI / 5. * point[0]);
        Bz(ix)     = std::sin(2 * M_PI / 5. * point[0]);
    }

    ampere.setLayout(&layout);
    ampere(B, J);

    auto psi_p = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    for (auto ix = psi_p; ix <= pei_p; ++ix)
    {
        EXPECT_THAT(Jy(ix), ::testing::DoubleNear((expectedJy[ix]), 1e-12));
        EXPECT_THAT(Jz(ix), ::testing::DoubleNear((expectedJz[ix]), 1e-12));
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
