

#ifndef PHARE_CORE_DEF_DETAIL_MKN_GPU_HPP
#define PHARE_CORE_DEF_DETAIL_MKN_GPU_HPP

#if !defined(PHARE_HAVE_MKN_GPU)
#define PHARE_HAVE_MKN_GPU 0
#endif

#if PHARE_HAVE_MKN_GPU
#include "mkn/gpu.hpp"
#define PHARE_WITH_MKN_GPU(...) __VA_ARGS__
#else
#define PHARE_WITH_MKN_GPU(...)
#endif

#endif /* PHARE_CORE_DEF_DETAIL_MKN_GPU_HPP */
