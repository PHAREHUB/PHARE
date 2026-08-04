
#include "core/def/phare_mpi.hpp"
#include "phare_core.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"

#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/coarse_cell_round_out.hpp"
#include "amr/data/field/refine/composite_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_composite_refiner.hpp"
#include "amr/data/field/refine/magnetic_refine_patch_strategy.hpp"
#include "amr/data/field/refine/adpt_magnetic_refine_patch_strategy.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"

#include "core/utilities/box/box.hpp"
#include "core/utilities/index/index.hpp"

#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/tbox/SAMRAIManager.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>



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
// case covers the fill-all-faces stage-1 kernel + the runtime tensor product of the 1-D rows.

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

std::string boxStr(SAMRAI::hier::Box const& box)
{
    auto const dim = box.getDim().getValue();
    std::string s  = "[(";
    for (unsigned short d = 0; d < dim; ++d)
        s += std::to_string(box.lower(d)) + (d + 1 < dim ? "," : "),(");
    for (unsigned short d = 0; d < dim; ++d)
        s += std::to_string(box.upper(d)) + (d + 1 < dim ? "," : ")]");
    return s;
}

// gtest's EXPECT_EQ cannot see PHARE::amr's operator== for SAMRAI boxes (it is not found by ADL
// from testing::internal), and a raw byte comparison would also compare BoxIds.
::testing::AssertionResult boxesEqual(SAMRAI::hier::Box const& lhs, SAMRAI::hier::Box const& rhs)
{
    if (PHARE::amr::operator==(lhs, rhs))
        return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure() << boxStr(lhs) << " != " << boxStr(rhs);
}
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


// B (Bx: primal-x normal, dual-y tangential): the two tangential (y) children of a shared
// (even-x) face are antisymmetric about the coarse value ⇒ sum = 2·C ⇒ ∇·B-neutral (stage-1
// property, unaffected by the magnetic round-out). Interior (odd-x, normal-direction) faces are
// now filled too — stage 1 of the ADPT prolongation fills ALL fine faces; ownership of divB-safety
// there moves to the stage-2 touch-up (not exercised by this kernel-only test). Assert the interior
// fill matches the same tensor-product stencil value the kernel uses (directionalInterp /
// directionalProlongation), computed independently of CompositeFieldRefiner's own orchestration.
TEST(magneticCompositeRefiner2D, sharedFaceTangentialChildrenAreDivBNeutral)
{
    using ImplYee2D = GridYee2D::implT;

    auto src_at = [](int ix, int iy) { return 1.0 * ix + 0.5 * iy + 0.1 * ix * iy; };
    std::array<QtyCentering, 2> centering{QtyCentering::primal, QtyCentering::dual}; // Bx

    auto rowX = [](int parity) {
        std::vector<WeightPoint<2>> row;
        if (parity == 0)
            row.push_back({Point<int, 2>{}, 1.0});
        else
            for (auto const& w :
                 ImplYee2D::directionalInterp<0, ImplYee2D::InterpDir::PrimalToDual, 2>())
                row.push_back(w);
        return row;
    };
    auto rowY = [](int parity) {
        std::vector<WeightPoint<2>> row;
        if (parity == 0)
            for (auto const& w : ImplYee2D::directionalProlongation<1, -1, 2>())
                row.push_back(w);
        else
            for (auto const& w : ImplYee2D::directionalProlongation<1, +1, 2>())
                row.push_back(w);
        return row;
    };
    auto expectedAt = [&](int Ix, int Iy, int px, int py) {
        double v = 0.0;
        for (auto const& wx : rowX(px))
            for (auto const& wy : rowY(py))
                v += wx.coef * wy.coef * src_at(Ix + wx.indexes[0], Iy + wy.indexes[1]);
        return v;
    };

    Grid2D src{"c", HybridQuantity::Scalar::Bx, 6u, 6u};
    Grid2D dst{"f", HybridQuantity::Scalar::Bx, std::array<std::uint32_t, 2>{12u, 12u}, NaN};
    for (std::size_t ix = 0; ix < src.shape()[0]; ++ix)
        for (std::size_t iy = 0; iy < src.shape()[1]; ++iy)
            src(ix, iy) = src_at(ix, iy);

    MagneticCompositeRefiner<GridYee2D, Grid2D, 2> refiner{};
    refiner.refineBox(src, dst, boxOf<2>({4, 4}, {5, 7}), centering, boxOf<2>({0, 0}, {11, 11}),
                      boxOf<2>({0, 0}, {5, 5}), ratio2(2));

    // fine x=4 is a shared (even) face, coarse anchor Ix=2
    for (int iy = 2; iy <= 3; ++iy)
    {
        double const c0 = dst(4, 2 * iy);
        double const c1 = dst(4, 2 * iy + 1);
        EXPECT_NEAR(c0 + c1, 2.0 * src_at(2, iy), 1e-12);
    }

    // fine x=5 is interior (odd-x, normal-direction), coarse anchor Ix=2: stage 1 now fills it
    // (no more single-owner skip). Assert filled (non-NaN) and matching the tensor-product
    // stencil value built from the primitive weight tables above.
    for (int iy = 4; iy <= 7; ++iy)
    {
        int const Iy = iy / 2;
        int const py = iy - 2 * Iy;
        EXPECT_FALSE(std::isnan(dst(5, iy)));
        EXPECT_NEAR(dst(5, iy), expectedAt(2, Iy, 1, py), 1e-12);
    }
}


// ----- the whole-coarse-cell round-out
// ------------------------------------------------------------

TEST(coarseCellRoundOut, parityHelpersAreFloorBasedOnNegativeIndices)
{
    // `i % 2` is -1 on odd negatives, so a `== 1` parity test would silently invert there. These
    // matter: the measured failing fill boxes really do have lower indices -3, -2, -1.
    EXPECT_TRUE(isOddIndex(-3));
    EXPECT_FALSE(isOddIndex(-2));
    EXPECT_TRUE(isOddIndex(-1));
    EXPECT_FALSE(isOddIndex(0));
    EXPECT_TRUE(isOddIndex(1));

    EXPECT_EQ(roundDownToEven(-3), -4);
    EXPECT_EQ(roundDownToEven(-2), -2);
    EXPECT_EQ(roundDownToEven(-1), -2);
    EXPECT_EQ(roundDownToEven(0), 0);
    EXPECT_EQ(roundDownToEven(1), 0);

    EXPECT_EQ(roundUpToOddIndex(-4), -3);
    EXPECT_EQ(roundUpToOddIndex(-3), -3);
    EXPECT_EQ(roundUpToOddIndex(-2), -1);
    EXPECT_EQ(roundUpToOddIndex(-1), -1);
    EXPECT_EQ(roundUpToOddIndex(0), 1);

    EXPECT_EQ(roundUpToEvenIndex(-3), -2);
    EXPECT_EQ(roundUpToEvenIndex(-2), -2);
    EXPECT_EQ(roundUpToEvenIndex(-1), 0);
    EXPECT_EQ(roundUpToEvenIndex(0), 0);
    EXPECT_EQ(roundUpToEvenIndex(1), 2);
}


TEST(coarseCellRoundOut, cellAndFieldBoxesAgreeOnWholeCoarseCells)
{
    // the measured failing fill box: one cell row, the lower half of coarse row 32
    auto const fill    = boxOf<2>({-2, 64}, {81, 64});
    auto const rounded = roundCellBoxOutToCoarseCells<2>(fill);
    EXPECT_TRUE(boxesEqual(rounded, boxOf<2>({-2, 64}, {81, 65})));

    // negative, odd on both lower edges (the measured coarse-interp temporary box)
    EXPECT_TRUE(boxesEqual(roundCellBoxOutToCoarseCells<2>(boxOf<2>({-1, 47}, {80, 63})),
                           boxOf<2>({-2, 46}, {81, 63})));

    // a field box rounds to the field box of the rounded cell box: dual rows 2I..2J+1, primal
    // nodes 2I..2J+2
    std::array<QtyCentering, 2> bxCentering{QtyCentering::primal, QtyCentering::dual}; // Bx
    EXPECT_TRUE(
        boxesEqual(roundFieldBoxOutToCoarseCells<2>(boxOf<2>({-2, 64}, {82, 64}), bxCentering),
                   boxOf<2>({-2, 64}, {82, 65})));

    std::array<QtyCentering, 2> byCentering{QtyCentering::dual, QtyCentering::primal}; // By
    EXPECT_TRUE(
        boxesEqual(roundFieldBoxOutToCoarseCells<2>(boxOf<2>({-2, 64}, {81, 65}), byCentering),
                   boxOf<2>({-2, 64}, {81, 66})));
}


TEST(coarseCellRoundOut, wholeCoarseCellsMeansLowerEvenAndUpperOdd)
{
    EXPECT_TRUE(isWholeCoarseCells<2>(boxOf<2>({-2, 64}, {81, 65})));
    EXPECT_TRUE(isWholeCoarseCells<2>(boxOf<2>({0, 0}, {1, 1})));      // a single coarse cell
    EXPECT_FALSE(isWholeCoarseCells<2>(boxOf<2>({-3, 64}, {81, 65}))); // lower x odd
    EXPECT_FALSE(isWholeCoarseCells<2>(boxOf<2>({-2, 64}, {80, 65}))); // upper x even
    EXPECT_FALSE(isWholeCoarseCells<2>(boxOf<2>({-2, 64}, {81, 64}))); // upper y even
    // a rounded-out box always satisfies it, negative indices included
    EXPECT_TRUE(isWholeCoarseCells<2>(roundCellBoxOutToCoarseCells<2>(boxOf<2>({-3, -1}, {4, 6}))));
}


// ----- magnetic prolongation over a misaligned fill box
// -------------------------------------------
//
// The regression these two guard. SAMRAI fill boxes are not unions of whole coarse cells — a
// recursive schedule's coarse-interpolation temporary is filled plus a one-cell ring, and one cell
// is always half a coarse cell. Tóth-Roe reaches exactly one coarse cell, so on such a box it used
// to read shared faces that lay outside the box and had therefore never been written: NaN B on the
// fine level, then NaN E, then a particle out of its patch.

namespace
{
// Coarse data that is *exactly* discretely divergence-free (unit spacing):
//   dBx/di + dBy/dj = A + (-A) = 0
// and non-trivial in both directions, so the Tóth-Roe transverse terms are exercised.
constexpr double coefA = 0.7, coefG = 0.35, coefD = -0.2, coefC1 = 1.1, coefC2 = -0.4;

double coarseBxAt(int I, int J)
{
    return coefA * I + coefG * J + coefC1;
}
double coarseByAt(int I, int J)
{
    return -coefA * J + coefD * I + coefC2;
}

// What order-2 prolongation followed by Tóth-Roe must produce. The dual ±¼ ladder is exact on a
// linear profile, so the shared face at (2I, 2J+p) is coarseBxAt(I,J) ± G/4; Tóth-Roe's base term
// is exact on a linear profile and its transverse term is a mixed second difference, which vanishes
// on data with no cross term. Both fine components are again linear, with fine divergence
// A/2 - A/2 = 0.
double fineBxAt(int i, int j)
{
    return 0.5 * coefA * i + 0.5 * coefG * j - 0.25 * coefG + coefC1;
}
double fineByAt(int i, int j)
{
    return -0.5 * coefA * j + 0.5 * coefD * i - 0.25 * coefD + coefC2;
}

constexpr int fieldGhosts = static_cast<int>(GridYee2D::options.field_ghost_width);

using VecFieldData2D = PHARE::amr::TensorFieldData<1, GridYee2D, Grid2D, HybridQuantity>;
using MagStrategy2D  = MagneticRefinePatchStrategy<int, VecFieldData2D>;
using FieldGeom2D    = FieldGeometry<GridYee2D, HybridQuantity::Scalar>;

GridYee2D layoutOf(SAMRAI::hier::Box const& cellBox, double dl)
{
    std::array<std::uint32_t, 2> nbrCells{static_cast<std::uint32_t>(cellBox.numberCells(0)),
                                          static_cast<std::uint32_t>(cellBox.numberCells(1))};
    return GridYee2D{{dl, dl}, nbrCells, Point<double, 2>{0., 0.}, phare_box_from<2>(cellBox)};
}

//! allocated field box of a layout, in AMR field indices
SAMRAI::hier::Box ghostFieldBoxOf(GridYee2D const& layout, HybridQuantity::Scalar qty)
{
    auto const cells = samrai_box_from(grow(layout.AMRBox(), fieldGhosts));
    return FieldGeom2D::toFieldBox(cells, qty, layout);
}

//! the local index of an AMR field index; identical for both centerings (local = amr - lower + g)
auto localOf(GridYee2D const& layout, int i, int j)
{
    return layout.AMRToLocal(Point<int, 2>{i, j});
}

std::array<Grid2D, 3> nanBFields(GridYee2D const& layout)
{
    return {Grid2D{"Bx", layout, HybridQuantity::Scalar::Bx, NaN},
            Grid2D{"By", layout, HybridQuantity::Scalar::By, NaN},
            Grid2D{"Bz", layout, HybridQuantity::Scalar::Bz, NaN}};
}
} // namespace


// A 1-cell-thick fill strip whose lower x is odd and whose single cell row is odd: every edge cuts
// a coarse cell in half. Round-out (gather side and postprocess side) must make every reconstructed
// interior face exist, be finite, and be div-free.
TEST(magneticProlongation2D, misalignedFillBoxReconstructsFiniteDivBFreeInteriorFaces)
{
    auto const finePatch   = boxOf<2>({0, 0}, {15, 15});
    auto const coarsePatch = boxOf<2>({0, 0}, {7, 7});

    auto const fineLayout   = layoutOf(finePatch, 0.1);
    auto const coarseLayout = layoutOf(coarsePatch, 0.2);

    // deliberately misaligned: lower x odd, upper x even, one cell row at odd y
    auto const fill = boxOf<2>({1, 3}, {12, 3});

    auto fields = nanBFields(fineLayout);

    // ---- gather the shared faces, exactly as KernelTensorFieldRefineOperator::refine does
    MagneticCompositeRefiner<GridYee2D, Grid2D, 2> refiner{};
    std::array bQties{HybridQuantity::Scalar::Bx, HybridQuantity::Scalar::By};
    for (std::size_t c = 0; c < bQties.size(); ++c)
    {
        Grid2D coarse{"c", coarseLayout, bQties[c], NaN};
        auto const srcFieldBox = ghostFieldBoxOf(coarseLayout, bQties[c]);
        for (std::uint32_t lx = 0; lx < coarse.shape()[0]; ++lx)
            for (std::uint32_t ly = 0; ly < coarse.shape()[1]; ++ly)
            {
                int const I    = static_cast<int>(lx) + srcFieldBox.lower(0);
                int const J    = static_cast<int>(ly) + srcFieldBox.lower(1);
                coarse(lx, ly) = c == 0 ? coarseBxAt(I, J) : coarseByAt(I, J);
            }

        auto const dstFieldBox = ghostFieldBoxOf(fineLayout, bQties[c]);
        auto const overlapBox  = FieldGeom2D::toFieldBox(fill, bQties[c], fineLayout);

        refiner.refineBox(coarse, fields[c], dstFieldBox * overlapBox,
                          fineLayout.centering(bQties[c]), dstFieldBox, srcFieldBox, ratio2(2));
    }

    // ---- reconstruct the interior faces, exactly as postprocessRefine does
    auto const region = MagStrategy2D::reconstructionRegion(
        fill, samrai_box_from(grow(fineLayout.AMRBox(), fieldGhosts)));
    EXPECT_TRUE(boxesEqual(region, boxOf<2>({0, 2}, {13, 3}))); // rounded out, no clip needed

    MagStrategy2D::reconstructInteriorFaces(fields, fineLayout, region);

    // every fine face of the region — shared (gathered) and interior (reconstructed) alike — is
    // finite and equal to the exact prolongation
    auto const& bx         = fields[0];
    auto const& by         = fields[1];
    std::size_t interiorBx = 0, interiorBy = 0;

    for (auto const& idx :
         phare_box_from<2>(FieldGeom2D::toFieldBox(region, HybridQuantity::Scalar::Bx, fineLayout)))
    {
        auto const loc = localOf(fineLayout, idx[0], idx[1]);
        ASSERT_TRUE(std::isfinite(bx(loc[0], loc[1])))
            << "Bx not written at amr (" << idx[0] << "," << idx[1] << ")";
        EXPECT_NEAR(bx(loc[0], loc[1]), fineBxAt(idx[0], idx[1]), 1e-12);
        if (idx[0] % 2 != 0)
            ++interiorBx;
    }

    for (auto const& idx :
         phare_box_from<2>(FieldGeom2D::toFieldBox(region, HybridQuantity::Scalar::By, fineLayout)))
    {
        auto const loc = localOf(fineLayout, idx[0], idx[1]);
        ASSERT_TRUE(std::isfinite(by(loc[0], loc[1])))
            << "By not written at amr (" << idx[0] << "," << idx[1] << ")";
        EXPECT_NEAR(by(loc[0], loc[1]), fineByAt(idx[0], idx[1]), 1e-12);
        if (idx[1] % 2 != 0)
            ++interiorBy;
    }

    EXPECT_GT(interiorBx, 0u); // the test is only meaningful if faces were reconstructed
    EXPECT_GT(interiorBy, 0u);

    // and the discrete divergence of every fine cell of the region is zero
    for (auto const& cell : phare_box_from<2>(region))
    {
        auto const bxLo = localOf(fineLayout, cell[0], cell[1]);
        auto const bxHi = localOf(fineLayout, cell[0] + 1, cell[1]);
        auto const byHi = localOf(fineLayout, cell[0], cell[1] + 1);

        double const divB = (bx(bxHi[0], bxHi[1]) - bx(bxLo[0], bxLo[1]))
                            + (by(byHi[0], byHi[1]) - by(bxLo[0], bxLo[1]));
        EXPECT_NEAR(divB, 0., 1e-12) << "divB at fine cell (" << cell[0] << "," << cell[1] << ")";
    }
}


// When the fill box already reaches the allocation, round-out is clipped back and can leave a
// coarse cell half-covered. There is no well-posed reconstruction for such a cell — some of its
// inputs are faces nothing ever wrote — so this is rejected rather than silently skipped. It cannot
// happen in any supported configuration (field_ghost_width >= fill_ring + 1 holds everywhere),
// which is why a deliberately misaligned *allocation* is needed to reach it at all.
TEST(magneticProlongation2D, halfCoveredCoarseCellsAreRejected)
{
    // a coarse-interpolation temporary: its box is coarsen(unfilled), so its parity is arbitrary.
    // lower x odd and upper x even ⇒ the allocated box is misaligned in x whatever the ghost width,
    // and a fill box reaching that edge cannot be rounded out inside the allocation.
    auto const fineLayout = layoutOf(boxOf<2>({-1, 0}, {14, 15}), 0.1);
    auto const ghostBox   = samrai_box_from(grow(fineLayout.AMRBox(), fieldGhosts));

    EXPECT_FALSE(isWholeCoarseCells<2>(ghostBox));
    EXPECT_THROW(MagStrategy2D::reconstructionRegion(ghostBox, ghostBox), std::runtime_error);

    // an aligned allocation of the same size is accepted, so it really is the parity that is
    // rejected
    auto const alignedLayout = layoutOf(boxOf<2>({0, 0}, {15, 15}), 0.1);
    auto const alignedGhosts = samrai_box_from(grow(alignedLayout.AMRBox(), fieldGhosts));
    EXPECT_NO_THROW(MagStrategy2D::reconstructionRegion(alignedGhosts, alignedGhosts));
}


// ==================================================================================================
// ADPTMagneticRefinePatchStrategy stage-2 touch-up tests (Balsara divB-free prolongation)
// ==================================================================================================
//
// Exercises the public static correctBx2d/correctBy2d directly: they are plain static functions,
// so no SAMRAI Patch/ResourcesManager machinery is needed. Stage 1 (CompositeFieldRefiner, reused
// from the value-level tests above) fills every fine face of Bx/By from coarse data; these statics
// then apply the stage-2 divergence-equalizing correction, sharing one DivCache per postprocess
// pass -- the same contract as ADPTMagneticRefinePatchStrategy::postprocessRefine.
//
// The strategy class only reads, from its TensorFieldDataT template parameter, the compile-time
// typedefs used at class scope (Geometry/gridlayout_type/N/dimension); Geometry and ResMan are
// never touched by the statics under test, so both are empty stand-ins.
//
// order-4 is deliberately NOT exercised here (unlike the upstream ADPT gtests this was ported
// from): CompositeFieldRefiner on this branch static_asserts order == 2 (higher-order-refinement
// caps refinement order at 2 for hybrid), so an order-4 fillFaces2D<4> would fail to compile. The
// touch-up statics themselves are order-independent, so order 2 alone still exercises the full
// stage-2 contract (memoised stage-1 snapshot, per-zone equalization).
namespace
{
struct DummyGeometry2D
{
};

struct DummyTensorFieldData2D
{
    using Geometry                         = DummyGeometry2D;
    using gridlayout_type                  = GridYee2D;
    static constexpr std::size_t N         = 3;
    static constexpr std::size_t dimension = 2;
};

struct DummyResMan2D
{
};

using ADPT2D = ADPTMagneticRefinePatchStrategy<DummyResMan2D, DummyTensorFieldData2D>;

// A GridLayout whose AMRToLocal is the identity (AMR index == array-local index): avoids
// depending on GridLayoutImplYee's internal ghost-width value, matching the "lower=0 =>
// AMR==local" convention the CompositeFieldRefiner value tests above already rely on.
GridYee2D identityLayout2D()
{
    using CoreBox = PHARE::core::Box<int, 2>;
    GridYee2D probe{{1., 1.}, {40u, 40u}, {{0., 0.}}, CoreBox{Point{0, 0}, Point{39, 39}}};
    auto const g0 = probe.AMRToLocal(Point{0, 0});
    int const G   = static_cast<int>(g0[0]);
    return GridYee2D{{1., 1.}, {40u, 40u}, {{0., 0.}}, CoreBox{Point{G, G}, Point{G + 39, G + 39}}};
}

// raw ("flux") divergence of fine cell (cx,cy): the same quantity the strategy's
// subzoneDiv2d_ computes internally.
double rawDiv2d(Grid2D& bx, Grid2D& by, int cx, int cy)
{
    return (bx(cx + 1, cy) - bx(cx, cy)) + (by(cx, cy + 1) - by(cx, cy));
}

Grid2D makeGrid2D(HybridQuantity::Scalar qty, int n)
{
    return Grid2D{"g", qty, static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(n)};
}

constexpr int coarseN_ = 12;
constexpr int fineN_   = 2 * coarseN_;

// fills every fine face of bx (primal-x/dual-y) and by (dual-x/primal-y) inside fine indices
// [6,17]^2 (coarse cells [3,8]) from full [0,coarseN_)^2 coarse arrays, via CompositeFieldRefiner.
template<int order>
void fillFaces2D(Grid2D& bxCoarse, Grid2D& byCoarse, Grid2D& bxFine, Grid2D& byFine)
{
    std::array<QtyCentering, 2> const bxCentering{QtyCentering::primal, QtyCentering::dual};
    std::array<QtyCentering, 2> const byCentering{QtyCentering::dual, QtyCentering::primal};

    auto const destBox   = boxOf<2>({6, 6}, {17, 17});
    auto const destGhost = boxOf<2>({0, 0}, {fineN_ - 1, fineN_ - 1});
    auto const srcGhost  = boxOf<2>({0, 0}, {coarseN_ - 1, coarseN_ - 1});

    CompositeFieldRefiner<GridYee2D, Grid2D, order> refiner{};
    refiner.refineBox(bxCoarse, bxFine, destBox, bxCentering, destGhost, srcGhost, ratio2(2));
    refiner.refineBox(byCoarse, byFine, destBox, byCentering, destGhost, srcGhost, ratio2(2));
}

// applies the stage-2 touch-up over the same fine box used to fill the faces (matches
// ADPTMagneticRefinePatchStrategy::postprocessRefine's per-component loop shape: the parity
// gate inside correctBx2d/correctBy2d selects only the interior/odd faces).
void touchUp2D(Grid2D& bxFine, Grid2D& byFine)
{
    auto const layout  = identityLayout2D();
    auto const destBox = boxOf<2>({6, 6}, {17, 17});

    ADPT2D::DivCache cache;
    for (auto const& i : phare_box_from<2>(destBox))
        ADPT2D::correctBx2d(cache, bxFine, byFine, layout, i);
    for (auto const& i : phare_box_from<2>(destBox))
        ADPT2D::correctBy2d(cache, bxFine, byFine, layout, i);
}

void fillNaN2D(Grid2D& bxFine, Grid2D& byFine)
{
    for (int i = 0; i < fineN_; ++i)
        for (int j = 0; j < fineN_; ++j)
        {
            bxFine(i, j) = NaN;
            byFine(i, j) = NaN;
        }
}
} // namespace


// Coarse B built from a stream function psi (bx = y-difference, by = -x-difference of psi at the
// same 4 corner nodes bounding the cell) is discretely divergence-free by construction, at ANY
// grid resolution, independently of psi's shape -- the corner-algebra telescopes to zero
// regardless of mesh spacing. psi here is a generic (non-affine) cubic, so the order-2 stage-1
// fill is not already exact: the stage-2 touch-up is what must zero the fine divergence.
template<int order>
static void runDivFreeCase()
{
    auto psi = [](int x, int y) { return 0.7 * x * x * y - 0.3 * x * y * y + 1.3 * x - 0.6 * y; };

    Grid2D bxCoarse = makeGrid2D(HybridQuantity::Scalar::Bx, coarseN_);
    Grid2D byCoarse = makeGrid2D(HybridQuantity::Scalar::By, coarseN_);
    Grid2D bxFine   = makeGrid2D(HybridQuantity::Scalar::Bx, fineN_);
    Grid2D byFine   = makeGrid2D(HybridQuantity::Scalar::By, fineN_);

    for (int I = 0; I < coarseN_; ++I)
        for (int J = 0; J < coarseN_; ++J)
        {
            bxCoarse(I, J) = psi(I, J + 1) - psi(I, J);
            byCoarse(I, J) = -(psi(I + 1, J) - psi(I, J));
        }
    fillNaN2D(bxFine, byFine);

    fillFaces2D<order>(bxCoarse, byCoarse, bxFine, byFine);
    touchUp2D(bxFine, byFine);

    for (int cx = 8; cx <= 15; ++cx)
        for (int cy = 8; cy <= 15; ++cy)
            EXPECT_NEAR(rawDiv2d(bxFine, byFine, cx, cy), 0.0, 1e-12)
                << "order=" << order << " cx=" << cx << " cy=" << cy;
}

TEST(ADPTMagneticTouchUp2D, correctsToExactDivBFreeOnGenericDivFreeCoarseData)
{
    runDivFreeCase<2>();
}


// Same flow but the coarse field is NOT divergence-free: the touch-up doesn't (and can't) zero
// out the divergence -- per the class docstring it *equalizes* the 4 fine-subzone divergences of
// the coarse cell they split (they all become the transported zone divergence q0 != 0). Assert
// equality within each complete coarse zone, not zero.
template<int order>
static void runEqualizeCase()
{
    auto bxOf = [](int I, int J) { return I * I - 0.3 * J + 0.2 * I * J; };
    auto byOf = [](int I, int J) { return 0.5 * J * J + 0.1 * I - I * J; };

    Grid2D bxCoarse = makeGrid2D(HybridQuantity::Scalar::Bx, coarseN_);
    Grid2D byCoarse = makeGrid2D(HybridQuantity::Scalar::By, coarseN_);
    Grid2D bxFine   = makeGrid2D(HybridQuantity::Scalar::Bx, fineN_);
    Grid2D byFine   = makeGrid2D(HybridQuantity::Scalar::By, fineN_);

    for (int I = 0; I < coarseN_; ++I)
        for (int J = 0; J < coarseN_; ++J)
        {
            bxCoarse(I, J) = bxOf(I, J);
            byCoarse(I, J) = byOf(I, J);
        }
    fillNaN2D(bxFine, byFine);

    fillFaces2D<order>(bxCoarse, byCoarse, bxFine, byFine);
    touchUp2D(bxFine, byFine);

    for (int I = 4; I <= 7; ++I)
        for (int J = 4; J <= 7; ++J)
        {
            double const d00 = rawDiv2d(bxFine, byFine, 2 * I, 2 * J);
            double const d10 = rawDiv2d(bxFine, byFine, 2 * I + 1, 2 * J);
            double const d01 = rawDiv2d(bxFine, byFine, 2 * I, 2 * J + 1);
            double const d11 = rawDiv2d(bxFine, byFine, 2 * I + 1, 2 * J + 1);

            EXPECT_NEAR(d10, d00, 1e-12) << "order=" << order << " I=" << I << " J=" << J;
            EXPECT_NEAR(d01, d00, 1e-12) << "order=" << order << " I=" << I << " J=" << J;
            EXPECT_NEAR(d11, d00, 1e-12) << "order=" << order << " I=" << I << " J=" << J;
        }
}

TEST(ADPTMagneticTouchUp2D, equalizesSubzoneDivergenceOnGenericNonDivFreeCoarseData)
{
    runEqualizeCase<2>();
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
