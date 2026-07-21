#ifndef PHARE_TEST_CORE_NUMERICS_BOUNDARY_CONDITION_HYBRID_BC_TEST_FIXTURES_HPP
#define PHARE_TEST_CORE_NUMERICS_BOUNDARY_CONDITION_HYBRID_BC_TEST_FIXTURES_HPP

#include "gtest/gtest.h"

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/patch_field_accessor.hpp"
#include "core/numerics/boundary_condition/boundary_condition_context.hpp"
#include "core/utilities/box/box.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"

namespace PHARE::core
{

// ─── Constants ───────────────────────────────────────────────────────────────

static constexpr std::uint32_t nCells      = 10u;
static constexpr std::uint32_t ghostWidth  = GridLayoutImplYee<1, 1>::ghost_width;
static constexpr double        interiorValue = 5.0;
static constexpr double        ghostSentinel = 99.0;

static constexpr std::uint32_t nCellsX2D = 10u;
static constexpr std::uint32_t nCellsY2D = 8u;

static constexpr std::uint32_t nCellsX3D = 10u;
static constexpr std::uint32_t nCellsY3D = 8u;
static constexpr std::uint32_t nCellsZ3D = 6u;


// ─── Type aliases ─────────────────────────────────────────────────────────────

using GridLayout1D = GridLayout<GridLayoutImplYee<1, 1>>;
using NdArray1D    = NdArrayVector<1, double>;
using Grid1D       = Grid<NdArray1D, HybridQuantity::Scalar>;
using Field1D      = Grid1D::field_type;
using VecField1D   = VecField<Field1D, HybridQuantity>;

using GridLayout2D = GridLayout<GridLayoutImplYee<2, 1>>;
using NdArray2D    = NdArrayVector<2, double>;
using Grid2D       = Grid<NdArray2D, HybridQuantity::Scalar>;
using Field2D      = Grid2D::field_type;
using VecField2D   = VecField<Field2D, HybridQuantity>;

using GridLayout3D = GridLayout<GridLayoutImplYee<3, 1>>;
using NdArray3D    = NdArrayVector<3, double>;
using Grid3D       = Grid<NdArray3D, HybridQuantity::Scalar>;
using Field3D      = Grid3D::field_type;
using VecField3D   = VecField<Field3D, HybridQuantity>;


// ─── Null accessor ────────────────────────────────────────────────────────────

template<typename FieldT>
struct NullFieldAccessorT : IPatchFieldAccessor<FieldT, HybridQuantity>
{
    FieldT& getField(HybridQuantity::Scalar) const override
    {
        throw std::runtime_error("NullFieldAccessorT: getField() should not be called");
    }
    VecField<FieldT, HybridQuantity> getVecField(HybridQuantity::Vector) const override
    {
        throw std::runtime_error("NullFieldAccessorT: getVecField() should not be called");
    }
};

using NullFieldAccessor = NullFieldAccessorT<Field1D>;


// Build a BC context bound to a single accessor; tests that don't care about the previous
// substage state pass the same accessor for both new and old. time=0 by default.
template<typename FieldT>
auto makeCtx(NullFieldAccessorT<FieldT> const& acc, double time = 0.0)
{
    return PHARE::core::BoundaryConditionContext<FieldT, HybridQuantity>{acc, acc, time};
}


// ─── 1D ghost-cell boxes ─────────────────────────────────────────────────────

inline Box<std::uint32_t, 1> lowerGhostCellBox()
{
    return {Point<std::uint32_t, 1>{0u}, Point<std::uint32_t, 1>{ghostWidth - 1u}};
}
inline Box<std::uint32_t, 1> upperGhostCellBox()
{
    return {Point<std::uint32_t, 1>{ghostWidth + nCells},
            Point<std::uint32_t, 1>{2u * ghostWidth + nCells - 1u}};
}


// ─── 2D ghost-cell boxes ─────────────────────────────────────────────────────

inline Box<std::uint32_t, 2> xLowerGhostCellBox2D()
{
    return {{0u, 0u}, {ghostWidth - 1u, nCellsY2D + 2u * ghostWidth - 1u}};
}
inline Box<std::uint32_t, 2> xUpperGhostCellBox2D()
{
    return {{ghostWidth + nCellsX2D, 0u},
            {2u * ghostWidth + nCellsX2D - 1u, nCellsY2D + 2u * ghostWidth - 1u}};
}
inline Box<std::uint32_t, 2> yLowerGhostCellBox2D()
{
    return {{0u, 0u}, {nCellsX2D + 2u * ghostWidth - 1u, ghostWidth - 1u}};
}
inline Box<std::uint32_t, 2> yUpperGhostCellBox2D()
{
    return {{0u, ghostWidth + nCellsY2D},
            {nCellsX2D + 2u * ghostWidth - 1u, 2u * ghostWidth + nCellsY2D - 1u}};
}


// ─── 3D ghost-cell boxes ─────────────────────────────────────────────────────

inline Box<std::uint32_t, 3> zLowerGhostCellBox3D()
{
    return {{0u, 0u, 0u},
            {nCellsX3D + 2u * ghostWidth - 1u, nCellsY3D + 2u * ghostWidth - 1u, ghostWidth - 1u}};
}
inline Box<std::uint32_t, 3> zUpperGhostCellBox3D()
{
    return {{0u, 0u, ghostWidth + nCellsZ3D},
            {nCellsX3D + 2u * ghostWidth - 1u, nCellsY3D + 2u * ghostWidth - 1u,
             2u * ghostWidth + nCellsZ3D - 1u}};
}


// ─── 1D scalar field fixture ──────────────────────────────────────────────────

struct FieldBC1D : testing::Test
{
    GridLayout1D layout{{0.1}, {nCells}, {0.0}};
    NullFieldAccessor acc;

    static constexpr auto qty = HybridQuantity::Scalar::rho;
    Grid1D grid{"rho", qty, layout.allocSize(qty)};
    Field1D& field{*(&grid)};

    std::uint32_t physStart{layout.physicalStartIndex(qty, Direction::X)};
    std::uint32_t physEnd{layout.physicalEndIndex(qty, Direction::X)};

    FieldBC1D()
    {
        for (std::uint32_t i = 0; i < grid.shape()[0]; ++i)
            field(i) = ghostSentinel;
        for (std::uint32_t i = physStart; i <= physEnd; ++i)
            field(i) = interiorValue;
    }
};


// ─── 1D VecField (B) fixture ──────────────────────────────────────────────────

struct VecFieldBC1D : testing::Test
{
    GridLayout1D layout{{0.1}, {nCells}, {0.0}};
    NullFieldAccessor acc;

    static constexpr auto vecQty = HybridQuantity::Vector::B;
    UsableTensorField<1, 1> B{"B", layout, vecQty};

    VecFieldBC1D()
    {
        for (std::size_t comp = 0; comp < 3; ++comp)
        {
            auto& f = B[comp];
            for (std::uint32_t i = 0; i < f.shape()[0]; ++i)
                f(i) = ghostSentinel;
            auto qty        = HybridQuantity::componentsQuantities(vecQty)[comp];
            std::uint32_t s = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t e = layout.physicalEndIndex(qty, Direction::X);
            for (std::uint32_t i = s; i <= e; ++i)
                f(i) = interiorValue;
        }
    }
};


// ─── 2D VecField (B) fixture ──────────────────────────────────────────────────

struct VecFieldBC2D : testing::Test
{
    GridLayout2D layout{{0.1, 0.1}, {nCellsX2D, nCellsY2D}, {0.0, 0.0}};
    NullFieldAccessorT<Field2D> acc;

    static constexpr auto vecQty = HybridQuantity::Vector::B;
    UsableTensorField<2, 1> B{"B", layout, vecQty};

    VecFieldBC2D()
    {
        for (std::size_t comp = 0; comp < 3; ++comp)
        {
            auto& f    = B[comp];
            auto shape = f.shape();
            for (std::uint32_t ix = 0; ix < shape[0]; ++ix)
                for (std::uint32_t iy = 0; iy < shape[1]; ++iy)
                    f(ix, iy) = ghostSentinel;
            auto qty         = HybridQuantity::componentsQuantities(vecQty)[comp];
            std::uint32_t sx = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t ex = layout.physicalEndIndex(qty, Direction::X);
            std::uint32_t sy = layout.physicalStartIndex(qty, Direction::Y);
            std::uint32_t ey = layout.physicalEndIndex(qty, Direction::Y);
            for (std::uint32_t ix = sx; ix <= ex; ++ix)
                for (std::uint32_t iy = sy; iy <= ey; ++iy)
                    f(ix, iy) = interiorValue;
        }
    }
};


// ─── 2D VecField (B) fixture with non-uniform By ─────────────────────────────

struct VecFieldBC2DNonUniformBy : testing::Test
{
    GridLayout2D layout{{0.1, 0.1}, {nCellsX2D, nCellsY2D}, {0.0, 0.0}};
    NullFieldAccessorT<Field2D> acc;

    static constexpr auto vecQty = HybridQuantity::Vector::B;
    UsableTensorField<2, 1> B{"B", layout, vecQty};

    VecFieldBC2DNonUniformBy()
    {
        for (std::size_t comp = 0; comp < 3; ++comp)
        {
            auto& f    = B[comp];
            auto shape = f.shape();
            for (std::uint32_t ix = 0; ix < shape[0]; ++ix)
                for (std::uint32_t iy = 0; iy < shape[1]; ++iy)
                    f(ix, iy) = ghostSentinel;

            auto qty         = HybridQuantity::componentsQuantities(vecQty)[comp];
            std::uint32_t sx = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t ex = layout.physicalEndIndex(qty, Direction::X);
            std::uint32_t sy = layout.physicalStartIndex(qty, Direction::Y);
            std::uint32_t ey = layout.physicalEndIndex(qty, Direction::Y);

            for (std::uint32_t ix = sx; ix <= ex; ++ix)
                for (std::uint32_t iy = sy; iy <= ey; ++iy)
                    f(ix, iy) = comp == 1 ? static_cast<double>(iy) : interiorValue;
        }
    }
};


// ─── 3D VecField (B) fixture ──────────────────────────────────────────────────

struct VecFieldBC3D : testing::Test
{
    GridLayout3D layout{{0.1, 0.1, 0.1}, {nCellsX3D, nCellsY3D, nCellsZ3D}, {0.0, 0.0, 0.0}};
    NullFieldAccessorT<Field3D> acc;

    static constexpr auto vecQty = HybridQuantity::Vector::B;
    UsableTensorField<3, 1> B{"B", layout, vecQty};

    VecFieldBC3D()
    {
        for (std::size_t comp = 0; comp < 3; ++comp)
        {
            auto& f    = B[comp];
            auto shape = f.shape();
            for (std::uint32_t ix = 0; ix < shape[0]; ++ix)
                for (std::uint32_t iy = 0; iy < shape[1]; ++iy)
                    for (std::uint32_t iz = 0; iz < shape[2]; ++iz)
                        f(ix, iy, iz) = ghostSentinel;
            auto qty         = HybridQuantity::componentsQuantities(vecQty)[comp];
            std::uint32_t sx = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t ex = layout.physicalEndIndex(qty, Direction::X);
            std::uint32_t sy = layout.physicalStartIndex(qty, Direction::Y);
            std::uint32_t ey = layout.physicalEndIndex(qty, Direction::Y);
            std::uint32_t sz = layout.physicalStartIndex(qty, Direction::Z);
            std::uint32_t ez = layout.physicalEndIndex(qty, Direction::Z);
            for (std::uint32_t ix = sx; ix <= ex; ++ix)
                for (std::uint32_t iy = sy; iy <= ey; ++iy)
                    for (std::uint32_t iz = sz; iz <= ez; ++iz)
                        f(ix, iy, iz) = interiorValue;
        }
    }
};

} // namespace PHARE::core

#endif // PHARE_TEST_CORE_NUMERICS_BOUNDARY_CONDITION_HYBRID_BC_TEST_FIXTURES_HPP
