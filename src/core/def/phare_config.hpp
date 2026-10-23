#ifndef PHARE_CORE_DEF_PHARE_CONFIG_HPP
#define PHARE_CORE_DEF_PHARE_CONFIG_HPP


#if !defined(PHARE_HAVE_GPU_HARDWARE)
#define PHARE_HAVE_GPU_HARDWARE 0
#endif

#include "core/def/detail/build_config.hpp"
#include "core/def/detail/umpire.hpp"
#include "core/def/detail/raja.hpp" // checks for umpire to know if gpu
#include "core/def/detail/mkn_avx.hpp"
#include "core/def/detail/mkn_gpu.hpp"
#include "core/def/kokkos_tools.hpp"

#include <cstdint>

#if !defined(PHARE_HAVE_GPU)
#define PHARE_HAVE_GPU PHARE_HAVE_GPU_HARDWARE
#endif // !defined(PHARE_HAVE_GPU)


#if PHARE_HAVE_GPU

#define _PHARE_DEV_FN_ __device__
#define _PHARE_HST_FN_ __host__
#define _PHARE_ALL_FN_ _PHARE_HST_FN_ _PHARE_DEV_FN_
#define PHARE_WITH_GPU(...) __VA_ARGS__

#else // !PHARE_HAVE_GPU

#define _PHARE_DEV_FN_
#define _PHARE_HST_FN_
#define _PHARE_ALL_FN_
#define PHARE_WITH_GPU(...)

#endif // PHARE_HAVE_GPU

#include "core/def/detail/thrust.hpp"

namespace PHARE
{
struct CompileOptions
{
    static constexpr bool WithUmpire   = PHARE_HAVE_UMPIRE;
    static constexpr bool WithMknGpu   = PHARE_HAVE_MKN_GPU;
    static constexpr bool WithMknGpuHW = PHARE_HAVE_MKN_GPU_HW;
    static constexpr bool WithMknAVX   = PHARE_HAVE_MKN_AVX;
    static constexpr bool WithRAJA     = PHARE_HAVE_RAJA;
    static constexpr bool WithThrust   = PHARE_HAVE_THRUST;
    static constexpr bool WithKokkos   = PHARE_HAVE_KOKKOS;
    static constexpr bool WithGpu      = PHARE_HAVE_GPU_HARDWARE;
};


enum class AllocatorMode : std::uint16_t {
    CPU = 0,
    GPU_UNIFIED,
    // GPU, // unified for now
};

} // namespace PHARE


#endif /*PHARE_CORE_DEF_PHARE_CONFIG_HPP*/
