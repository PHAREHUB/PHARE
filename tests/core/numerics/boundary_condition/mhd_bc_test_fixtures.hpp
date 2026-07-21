#ifndef PHARE_TEST_CORE_NUMERICS_BOUNDARY_CONDITION_MHD_BC_TEST_FIXTURES_HPP
#define PHARE_TEST_CORE_NUMERICS_BOUNDARY_CONDITION_MHD_BC_TEST_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/patch_field_accessor.hpp"
#include "core/numerics/boundary_condition/boundary_condition_context.hpp"
#include "core/utilities/box/box.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

namespace PHARE::core
{

// MHD Yee layout with reconstruction_nghosts=1:
//   ghost_width = roundUpToEven(1+2) = 4  (same for all spatial dims)
static constexpr std::uint32_t mhdGhostWidth = GridLayoutImplYeeMHD<1, 1, 1>::ghost_width;

/**
 * @brief Concrete patch-field accessor for MHD unit tests (dim-templated).
 *
 * Holds references to the scalar grids and vector field fixtures that the
 * FieldTotalEnergyFromPressureBoundaryCondition needs at runtime.
 * getVecField() copy-constructs a VecField view from the underlying
 * UsableVecFieldMHD; buffer pointers are preserved so writes through the
 * returned view affect the actual data.
 */
template<std::size_t dim>
struct MHDPatchFieldAccessorTest : IPatchFieldAccessor<FieldMHD<dim>, MHDQuantity>
{
    using GridMHD        = Grid<NdArrayVector<dim, double>, MHDQuantity::Scalar>;
    using UsableVecField = UsableVecFieldMHD<dim>;
    using VecField_      = VecFieldMHD<dim>;

    GridMHD& rho;
    GridMHD& P;
    GridMHD& Etot;
    UsableVecField& rhoV;
    UsableVecField& Bvec;

    MHDPatchFieldAccessorTest(GridMHD& rho_, GridMHD& P_, GridMHD& Etot_, UsableVecField& rhoV_,
                              UsableVecField& Bvec_)
        : rho{rho_}
        , P{P_}
        , Etot{Etot_}
        , rhoV{rhoV_}
        , Bvec{Bvec_}
    {
    }

    FieldMHD<dim>& getField(MHDQuantity::Scalar qty) const override
    {
        switch (qty)
        {
            case MHDQuantity::Scalar::rho: return *(&rho);
            case MHDQuantity::Scalar::P: return *(&P);
            case MHDQuantity::Scalar::Etot: return *(&Etot);
            default: throw std::runtime_error("MHDPatchFieldAccessorTest: unsupported scalar qty");
        }
    }

    VecField_ getVecField(MHDQuantity::Vector qty) const override
    {
        switch (qty)
        {
            case MHDQuantity::Vector::rhoV: return rhoV.super();
            case MHDQuantity::Vector::B: return Bvec.super();
            default: throw std::runtime_error("MHDPatchFieldAccessorTest: unsupported vector qty");
        }
    }
};


// Build a BC context bound to a single accessor; tests that don't care about the previous
// substage state pass the same accessor for both new and old. time=0 by default.
template<std::size_t dim>
auto makeCtx(MHDPatchFieldAccessorTest<dim> const& acc, double time = 0.0)
{
    return BoundaryConditionContext<FieldMHD<dim>, MHDQuantity>{acc, acc, time};
}


// ─── 1D MHD types and ghost boxes ────────────────────────────────────────────

using GridLayoutMHD1D = GridLayout<GridLayoutImplYeeMHD<1, 1, 1>>;
using GridMHD1D       = Grid<NdArrayVector<1, double>, MHDQuantity::Scalar>;

static constexpr std::uint32_t nCellsMHD = 10u;

inline Box<std::uint32_t, 1> mhdLowerGhostCellBox()
{
    return {Point<std::uint32_t, 1>{0u}, Point<std::uint32_t, 1>{mhdGhostWidth - 1u}};
}
inline Box<std::uint32_t, 1> mhdUpperGhostCellBox()
{
    return {Point<std::uint32_t, 1>{mhdGhostWidth + nCellsMHD},
            Point<std::uint32_t, 1>{2u * mhdGhostWidth + nCellsMHD - 1u}};
}


// ─── 2D MHD types and ghost boxes ────────────────────────────────────────────

using GridLayoutMHD2D = GridLayout<GridLayoutImplYeeMHD<2, 1, 1>>;
using GridMHD2D       = Grid<NdArrayVector<2, double>, MHDQuantity::Scalar>;

static constexpr std::uint32_t nCellsMHDX2D = 10u;
static constexpr std::uint32_t nCellsMHDY2D = 8u;

inline Box<std::uint32_t, 2> mhd2DXLowerGhostBox()
{
    return {{0u, 0u}, {mhdGhostWidth - 1u, nCellsMHDY2D + 2u * mhdGhostWidth - 1u}};
}
inline Box<std::uint32_t, 2> mhd2DXUpperGhostBox()
{
    return {{mhdGhostWidth + nCellsMHDX2D, 0u},
            {2u * mhdGhostWidth + nCellsMHDX2D - 1u, nCellsMHDY2D + 2u * mhdGhostWidth - 1u}};
}
inline Box<std::uint32_t, 2> mhd2DYLowerGhostBox()
{
    return {{0u, 0u}, {nCellsMHDX2D + 2u * mhdGhostWidth - 1u, mhdGhostWidth - 1u}};
}
inline Box<std::uint32_t, 2> mhd2DYUpperGhostBox()
{
    return {{0u, mhdGhostWidth + nCellsMHDY2D},
            {nCellsMHDX2D + 2u * mhdGhostWidth - 1u, 2u * mhdGhostWidth + nCellsMHDY2D - 1u}};
}


// ─── 3D MHD types and ghost boxes ────────────────────────────────────────────

using GridLayoutMHD3D = GridLayout<GridLayoutImplYeeMHD<3, 1, 1>>;
using GridMHD3D       = Grid<NdArrayVector<3, double>, MHDQuantity::Scalar>;

static constexpr std::uint32_t nCellsMHDX3D = 10u;
static constexpr std::uint32_t nCellsMHDY3D = 8u;
static constexpr std::uint32_t nCellsMHDZ3D = 6u;

inline Box<std::uint32_t, 3> mhd3DXLowerGhostBox()
{
    return {{0u, 0u, 0u},
            {mhdGhostWidth - 1u, nCellsMHDY3D + 2u * mhdGhostWidth - 1u,
             nCellsMHDZ3D + 2u * mhdGhostWidth - 1u}};
}
inline Box<std::uint32_t, 3> mhd3DXUpperGhostBox()
{
    return {{mhdGhostWidth + nCellsMHDX3D, 0u, 0u},
            {2u * mhdGhostWidth + nCellsMHDX3D - 1u, nCellsMHDY3D + 2u * mhdGhostWidth - 1u,
             nCellsMHDZ3D + 2u * mhdGhostWidth - 1u}};
}
inline Box<std::uint32_t, 3> mhd3DYLowerGhostBox()
{
    return {{0u, 0u, 0u},
            {nCellsMHDX3D + 2u * mhdGhostWidth - 1u, mhdGhostWidth - 1u,
             nCellsMHDZ3D + 2u * mhdGhostWidth - 1u}};
}
inline Box<std::uint32_t, 3> mhd3DYUpperGhostBox()
{
    return {{0u, mhdGhostWidth + nCellsMHDY3D, 0u},
            {nCellsMHDX3D + 2u * mhdGhostWidth - 1u, 2u * mhdGhostWidth + nCellsMHDY3D - 1u,
             nCellsMHDZ3D + 2u * mhdGhostWidth - 1u}};
}
inline Box<std::uint32_t, 3> mhd3DZLowerGhostBox()
{
    return {{0u, 0u, 0u},
            {nCellsMHDX3D + 2u * mhdGhostWidth - 1u, nCellsMHDY3D + 2u * mhdGhostWidth - 1u,
             mhdGhostWidth - 1u}};
}
inline Box<std::uint32_t, 3> mhd3DZUpperGhostBox()
{
    return {{0u, 0u, mhdGhostWidth + nCellsMHDZ3D},
            {nCellsMHDX3D + 2u * mhdGhostWidth - 1u, nCellsMHDY3D + 2u * mhdGhostWidth - 1u,
             2u * mhdGhostWidth + nCellsMHDZ3D - 1u}};
}

} // namespace PHARE::core

#endif // PHARE_TEST_CORE_NUMERICS_BOUNDARY_CONDITION_MHD_BC_TEST_FIXTURES_HPP
