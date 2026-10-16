#ifndef PHARE_CORE_DEF_DETAIL_MKN_AVX_HPP
#define PHARE_CORE_DEF_DETAIL_MKN_AVX_HPP

#if !defined(PHARE_HAVE_MKN_AVX)
#if __has_include("mkn/avx.hpp")
#define PHARE_HAVE_MKN_AVX 1

#else // __has_include("mkn/avx.hpp")
#define PHARE_HAVE_MKN_AVX 0

#endif // __has_include("mkn/avx.hpp")
#endif // !defined(PHARE_HAVE_MKN_AVX)

#if PHARE_HAVE_MKN_AVX

#include "mkn/avx.hpp"
#define PHARE_WITH_MKN_AVX(...) __VA_ARGS__

#else

#define PHARE_WITH_MKN_AVX(...)

#endif // PHARE_HAVE_MKN_AVX

#endif /* PHARE_CORE_DEF_DETAIL_MKN_AVX_HPP */
