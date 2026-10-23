#ifndef PHARE_CORE_DEF_DETAIL_MKN_GPU_HPP
#define PHARE_CORE_DEF_DETAIL_MKN_GPU_HPP


// #if !defined(PHARE_HAVE_MKN_GPU) && __has_include("mkn/gpu.hpp")
// #undef PHARE_HAVE_MKN_GPU
// #endif


#if !defined(PHARE_HAVE_MKN_GPU) || PHARE_HAVE_MKN_GPU == 1
#if __has_include("mkn/gpu.hpp")
#include "mkn/gpu.hpp"
#include "mkn/kul/alloc.hpp"
#define PHARE_HAVE_MKN_GPU 1

#define PHARE_HAVE_MKN_GPU_HW (MKN_GPU_ROCM || MKN_GPU_CUDA)
#endif // __has_include("mkn/gpu.hpp")
#endif // !defined(PHARE_HAVE_MKN_GPU)


#if !defined(PHARE_HAVE_MKN_GPU)
#define PHARE_HAVE_MKN_GPU 0
#endif // PHARE_HAVE_MKN_GPU

#if !defined(PHARE_HAVE_MKN_GPU_HW)
#define PHARE_HAVE_MKN_GPU_HW 0
#endif // PHARE_HAVE_MKN_GPU_HW


#if PHARE_HAVE_MKN_GPU
#define PHARE_WITH_MKN_GPU(...) __VA_ARGS__
#endif // PHARE_HAVE_MKN_GPU

#if !defined(PHARE_WITH_MKN_GPU)
#define PHARE_WITH_MKN_GPU(...)
#endif

#endif /* PHARE_CORE_DEF_DETAIL_MKN_GPU_HPP */
