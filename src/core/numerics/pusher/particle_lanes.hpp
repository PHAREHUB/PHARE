#ifndef PHARE_CORE_PUSHER_PARTICLE_LANES_HPP
#define PHARE_CORE_PUSHER_PARTICLE_LANES_HPP

// Lane blocks of the hybrid push and deposit kernels: the ParticleArray stays an array of
// structures, the kernels stage blocks of laneCount particles into lanes ([component][lane]
// arrays or SIMD registers) and run the same arithmetic on all lanes.
//
// Backend (compile time, hyb_xsimd.hpp PHARE_HYB_XSIMD): explicit xsimd batches of 8 doubles
// (AVX-512), or the AutoVec lane loops left to the compiler.
// Push: blocks of 8 (boris.hpp pushBlock_); the xsimd field gather loads one contiguous row of
// 4 doubles per lane and stencil row, transposed in registers (hyb_xsimd.hpp loadRows).
// Deposit: per particle at interpolation order 1, by blocks (start indices and weights in lane
// loops, then the lanes deposited in particle order) from order 2 on.
// Cell kernels (xsimd backend, domain arrays kept sorted by cell, IonUpdater particle_sort):
// blocks are formed per cell row; a block whose lanes all stay in the row's cell at mid step
// interpolates the fields from the cell's stencil union, loaded once per cell (boris.hpp
// pushBlockCell_), and the deposit accumulates in a lane-private window per cell
// (interpolator/cell_deposit.hpp). Below cellMinPpc particles per cell the blocks of a cell are
// mostly partial, and the kernels above walk the rows instead (IonUpdater::cellKernels_).

#include <algorithm>
#include <array>
#include <cstddef>


namespace PHARE::core::hyb
{

// cell kernels for a sorted domain array only from this mean number of particles per cell
// (IonUpdater::cellKernels_), the choice made again every cellChoiceEvery correctors
inline constexpr int cellMinPpc           = 12;
inline constexpr unsigned cellChoiceEvery = 20;

// blocked deposit from this interpolation order on (order 1: per particle, faster for its
// 2^dim-node stencil)
template<std::size_t interpOrder>
inline constexpr bool blockedDeposit = interpOrder >= 2;

// vectorisation backend of the lane kernels (template parameter): the lane loops left to the
// compiler, or explicit xsimd batches (hyb_xsimd.hpp, which sets hyb::Backend)
struct AutoVec
{
};
struct Xsimd
{
};


inline constexpr std::size_t laneCount = 8;

// lane pointers of one block of laneCount particles: same member names as the push kernel's
// lane block (delta, cell, v), each pointing to laneCount contiguous values
template<std::size_t dim>
struct LanePtrs
{
    double* delta[dim];
    int* cell[dim];
    double* v[3];
    double* weight; // read by the deposit only
};


/** one block of laneCount particles, component-major: the lanes of the blocked and cell
 * deposits */
template<std::size_t dim>
struct alignas(64) LaneBlock
{
    double weight[laneCount];
    int cell[dim][laneCount];
    double delta[dim][laneCount];
    double v[3][laneCount];
};


/** AoS particles seen as a deposit source: the lanes of a block are staged into a caller's
 * LaneBlock (last particle repeated in unused lanes)
 */
template<std::size_t dim, typename Particles>
struct AoSLanes
{
    using Block = LaneBlock<dim>;

    Particles const& particles;
    std::size_t size;

    Block const& depositLanes(std::size_t const i0, Block& blk) const
    {
        for (std::size_t l = 0; l < laneCount; ++l)
        {
            auto const& part = particles[std::min(i0 + l, size - 1)];
            blk.weight[l]    = part.weight;
            for (std::size_t d = 0; d < dim; ++d)
            {
                blk.cell[d][l]  = part.iCell[d];
                blk.delta[d][l] = part.delta[d];
            }
            for (std::size_t c = 0; c < 3; ++c)
                blk.v[c][l] = part.v[c];
        }
        return blk;
    }
};


} // namespace PHARE::core::hyb

#endif
