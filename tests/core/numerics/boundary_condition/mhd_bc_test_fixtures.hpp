#ifndef PHARE_TEST_CORE_NUMERICS_BOUNDARY_CONDITION_MHD_BC_TEST_FIXTURES_HPP
#define PHARE_TEST_CORE_NUMERICS_BOUNDARY_CONDITION_MHD_BC_TEST_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/models/options/mhd_options.hpp"
#include "phare_simulator_options.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/utilities/box/box.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

namespace PHARE::core
{

template<std::size_t dim>
inline constexpr PHARE::SimOpts mhdSimOptions{
    .dimension           = dim,
    .interp_order        = 1,
    .reconstruction_type = PHARE::MHDOpts::ReconstructionType::Constant,
    .slope_limiter_type  = PHARE::MHDOpts::SlopeLimiterType::None,
    .riemann_solver_type = PHARE::MHDOpts::RiemannSolverType::Rusanov};

template<std::size_t dim>
inline constexpr PHARE::MHDFieldOptions<mhdSimOptions<dim>> mhdFieldOptions{};

static constexpr std::uint32_t mhdGhostWidth = mhdFieldOptions<1>.field_ghost_width;

using GridLayoutMHD1D                    = GridLayout<PHARE::MHDOptions<mhdFieldOptions<1>>{}>;
using GridMHD1D                          = Grid<NdArrayVector<1, double>, MHDQuantity::Scalar>;
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
using GridLayoutMHD2D                       = GridLayout<PHARE::MHDOptions<mhdFieldOptions<2>>{}>;
using GridMHD2D                             = Grid<NdArrayVector<2, double>, MHDQuantity::Scalar>;
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
using GridLayoutMHD3D                       = GridLayout<PHARE::MHDOptions<mhdFieldOptions<3>>{}>;
using GridMHD3D                             = Grid<NdArrayVector<3, double>, MHDQuantity::Scalar>;
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
#endif
