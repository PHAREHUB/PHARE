
#include "core/def/phare_mpi.hpp"
#include "phare_core.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"

#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/composite_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_composite_refiner.hpp"

#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/tbox/SAMRAIManager.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <limits>



using namespace PHARE::core;
using namespace PHARE::amr;

using testing::Eq;



template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct aFieldRefineOperator : public ::testing::Test
{
};

using aFieldRefineOperatorInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<1>, InterpConst<2>>,
                     std::pair<DimConst<1>, InterpConst<3>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<2>, InterpConst<2>>, std::pair<DimConst<2>, InterpConst<3>>,
                     std::pair<DimConst<3>, InterpConst<1>>, std::pair<DimConst<3>, InterpConst<2>>,
                     std::pair<DimConst<3>, InterpConst<3>>>;

TYPED_TEST_SUITE(aFieldRefineOperator, aFieldRefineOperatorInfos);


// instantiation gate: forces full compilation of CompositeFieldRefiner<...,order> (vtable →
// refineBox) and the additive KernelFieldRefineOperator across all dim/interp, order 2.
TYPED_TEST(aFieldRefineOperator, kernelRefineOperatorCanBeCreated)
{
    static constexpr auto dim    = typename TypeParam::first_type{}();
    static constexpr auto interp = typename TypeParam::second_type{}();

    using GridYee = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;
    using GridT   = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;

    auto linearKernel = makeRefineKernel<GridYee, GridT>(2);
    EXPECT_NE(linearKernel, nullptr);

    KernelFieldRefineOperator<GridYee, GridT> kernelRefine{linearKernel};

    auto magLinearKernel = makeMagneticRefineKernel<GridYee, GridT>(2);
    EXPECT_NE(magLinearKernel, nullptr);

    EXPECT_ANY_THROW((makeRefineKernel<GridYee, GridT>(3)));
    EXPECT_ANY_THROW((makeRefineKernel<GridYee, GridT>(4)));
    EXPECT_ANY_THROW((makeMagneticRefineKernel<GridYee, GridT>(0)));
    EXPECT_ANY_THROW((makeMagneticRefineKernel<GridYee, GridT>(4)));
}




// ----- value-level refineBox tests for the composite kernel --------------------------------------
//
// Boxes are placed with lower=0 so AMR == local indexing; ratio 2. The fine destination is filled
// with NaN so the kernel's NaN-guard writes every targeted index. The numeric core is separable, so
// 1-D exercises both primitives (dual ±¼ ladder, primal half-point); one 2-D magnetic
// case covers the sharedFacesOnly skip + the runtime tensor product of the 1-D rows.

namespace
{
using GridYee1D = PHARE::core::PHARE_Types<PHARE::SimOpts{1, 1}>::Hybrid::GridLayout_t;
using Grid1D    = Grid<NdArrayVector<1>, HybridQuantity::Scalar>;
using GridYee2D = PHARE::core::PHARE_Types<PHARE::SimOpts{2, 1}>::Hybrid::GridLayout_t;
using Grid2D    = Grid<NdArrayVector<2>, HybridQuantity::Scalar>;

template<std::size_t dim>
SAMRAI::hier::Box boxOf(std::array<int, dim> lo, std::array<int, dim> up)
{
    SAMRAI::tbox::Dimension d{dim};
    SAMRAI::hier::Index loi{d, 0}, upi{d, 0};
    for (std::size_t k = 0; k < dim; ++k)
    {
        loi(k) = lo[k];
        upi(k) = up[k];
    }
    return SAMRAI::hier::Box{loi, upi, SAMRAI::hier::BlockId{0}};
}

SAMRAI::hier::IntVector ratio2(std::size_t dim)
{
    return SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{static_cast<unsigned short>(dim)}, 2};
}

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
} // namespace


// the two fine children of a coarse dual cell always mean back to its average (conservation),
// at every order.
TEST(compositeRefiner1D, dualChildrenConserveCoarseAverage)
{
    std::array<double, 8> coarse = {1.0, 1.3, 0.4, 2.0, -0.5, 0.7, 0.9, 0.2};
    std::array<QtyCentering, 1> centering{QtyCentering::dual};

    auto run = [&](auto refiner) {
        Grid1D src{"c", HybridQuantity::Scalar::rho, 8u};
        Grid1D dst{"f", HybridQuantity::Scalar::rho, std::array<std::uint32_t, 1>{16u}, NaN};
        for (std::size_t i = 0; i < src.shape()[0]; ++i)
            src(i) = coarse[i];

        refiner.refineBox(src, dst, boxOf<1>({4}, {11}), centering, boxOf<1>({0}, {15}),
                          boxOf<1>({0}, {7}), ratio2(1));

        for (int ix = 2; ix <= 5; ++ix)
            EXPECT_NEAR(0.5 * (dst(2 * ix) + dst(2 * ix + 1)), coarse[ix], 1e-12);
    };

    run(CompositeFieldRefiner<GridYee1D, Grid1D, 2>{});
}


// row sum ≡ 1 (partition of unity, S7): on a constant coarse field every fine child equals the
// constant.
TEST(compositeRefiner1D, dualRowSumsToOne)
{
    constexpr double konst = 3.7;
    std::array<QtyCentering, 1> centering{QtyCentering::dual};

    auto run = [&](auto refiner) {
        Grid1D src{"c", HybridQuantity::Scalar::rho, 8u};
        Grid1D dst{"f", HybridQuantity::Scalar::rho, std::array<std::uint32_t, 1>{16u}, NaN};
        for (std::size_t i = 0; i < src.shape()[0]; ++i)
            src(i) = konst;
        refiner.refineBox(src, dst, boxOf<1>({4}, {11}), centering, boxOf<1>({0}, {15}),
                          boxOf<1>({0}, {7}), ratio2(1));
        for (int f = 4; f <= 11; ++f)
            EXPECT_NEAR(dst(f), konst, 1e-14);
    };

    run(CompositeFieldRefiner<GridYee1D, Grid1D, 2>{});
}


// B (Bx: primal-x normal, dual-y tangential), sharedFacesOnly: the two tangential (y) children of a
// shared (even-x) face are antisymmetric about the coarse value ⇒ sum = 2·C ⇒ ∇·B-neutral.
// Interior (odd-x) faces are left untouched (NaN) by the Tóth-Roe owner.
TEST(magneticCompositeRefiner2D, sharedFaceTangentialChildrenAreDivBNeutral)
{
    auto src_at = [](int ix, int iy) { return 1.0 * ix + 0.5 * iy + 0.1 * ix * iy; };
    std::array<QtyCentering, 2> centering{QtyCentering::primal, QtyCentering::dual}; // Bx

    auto run = [&](auto refiner) {
        Grid2D src{"c", HybridQuantity::Scalar::Bx, 6u, 6u};
        Grid2D dst{"f", HybridQuantity::Scalar::Bx, std::array<std::uint32_t, 2>{12u, 12u}, NaN};
        for (std::size_t ix = 0; ix < src.shape()[0]; ++ix)
            for (std::size_t iy = 0; iy < src.shape()[1]; ++iy)
                src(ix, iy) = src_at(ix, iy);

        refiner.refineBox(src, dst, boxOf<2>({4, 4}, {5, 7}), centering, boxOf<2>({0, 0}, {11, 11}),
                          boxOf<2>({0, 0}, {5, 5}), ratio2(2));

        // fine x=4 is a shared (even) face, coarse anchor ix=2; x=5 is interior (odd) → skipped
        for (int iy = 2; iy <= 3; ++iy)
        {
            double const c0 = dst(4, 2 * iy);
            double const c1 = dst(4, 2 * iy + 1);
            EXPECT_NEAR(c0 + c1, 2.0 * src_at(2, iy), 1e-12);
        }
        EXPECT_TRUE(std::isnan(dst(5, 4))); // interior-normal face untouched
    };

    run(MagneticCompositeRefiner<GridYee2D, Grid2D, 2>{});
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();


    int testResult = RUN_ALL_TESTS();

    // Finalize
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
