#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include <cassert>
#include <stdexcept>
#include <string_view>

#if defined(PHARE_WITH_GPU) && PHARE_WITH_GPU

#define _PHARE_DEV_FN_ __device__
#define _PHARE_HST_FN_ __host__
#define _PHARE_ALL_FN_ _PHARE_HST_FN_ _PHARE_DEV_FN_

namespace PHARE
{
// !WARNING!
// PGIC++/NVC++ CANNOT HANDLE separate device/host functions!
template<typename T>
inline void throw_runtime_error(T const&) __device__
{
    // gpu cannot throw
    assert(false);
}

template<typename T>
inline void throw_runtime_error(T const& err) __host__
{
    throw std::runtime_error(err);
}
} // namespace PHARE

#else

#define _PHARE_DEV_FN_
#define _PHARE_HST_FN_
#define _PHARE_ALL_FN_

namespace PHARE
{
template<typename T>
inline void throw_runtime_error(T const& err)
{
    throw std::runtime_error(err);
}
} // namespace PHARE

#endif

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO)
#define PHARE_DEBUG_DO(...) __VA_ARGS__
#else
#define PHARE_DEBUG_DO(...)
#endif


#define _PHARE_TO_STR(x) #x // convert macro text to string
#define PHARE_TO_STR(x) _PHARE_TO_STR(x)



#endif /* PHARE_CORE_DEF_HPP */
