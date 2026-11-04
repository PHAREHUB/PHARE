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

#include <array>
#include <cstddef>


namespace PHARE::core::hyb
{

// vectorisation backend of the lane kernels (template parameter): the lane loops left to the
// compiler, or explicit xsimd batches (hyb_xsimd.hpp, which sets hyb::Backend)
struct AutoVec
{
};
struct Xsimd
{
};


inline constexpr std::size_t laneCount = 8;


} // namespace PHARE::core::hyb

#endif
