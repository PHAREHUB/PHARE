#ifndef PHARE_CORE_OPERATORS_HPP
#define PHARE_CORE_OPERATORS_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"

#ifndef PHARE_HAVE_GPU
#define PHARE_HAVE_GPU 0
#endif

#include <atomic>

// device atomics (atomicAdd...) only exist in the device compilation pass; host calls of
// _PHARE_ALL_FN_ code (e.g. GPU_UNIFIED CPU fallbacks) use the std::atomic paths instead
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
#define PHARE_DEVICE_PASS 1
#else
#define PHARE_DEVICE_PASS 0
#endif


namespace PHARE::core
{
template<typename T, bool atomic, bool GPU>
struct Operators
{
    T static constexpr ONE = 1;

    static constexpr bool device_atomic = GPU and atomic and PHARE_DEVICE_PASS;

    static_assert(not std::is_const_v<T>); // doesn't make sense

    void operator+=(T const& v) _PHARE_ALL_FN_
    {
        if constexpr (device_atomic)
        {
            atomicAdd(&t, v);
        }
        else if constexpr (atomic)
        {
            auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(&t);
            T tmp          = atomic_t.load();
            while (!atomic_t.compare_exchange_weak(tmp, tmp + v)) {}
        }
        else
            t += v;
    }
    void operator+=(T const&& v) _PHARE_ALL_FN_ { (*this) += v; }

    void operator-=(T const& v) _PHARE_ALL_FN_
    {
        if constexpr (device_atomic)
        {
            atomicSub(&t, v);
        }
        else if constexpr (atomic)
        {
            auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(&t);
            T tmp          = atomic_t.load();
            while (!atomic_t.compare_exchange_weak(tmp, tmp - v)) {}
        }
        else
            t -= v;
    }
    void operator-=(T const&& v) _PHARE_ALL_FN_ { (*this) -= v; }

    auto increment_return_old() _PHARE_ALL_FN_ // postfix increment
    {
        if constexpr (device_atomic)
        {
            auto o = atomicAdd(&t, ONE);
            PHARE_ASSERT(o < t);
            return o;
        }
        else if constexpr (atomic)
        {
            return std::atomic_ref<T>{t}.fetch_add(ONE);
        }
        else
        {
            T tmp = t;
            ++t;
            return tmp;
        }
    }

    auto static compare_and_swap(T* addr, T compare, T value) _PHARE_ALL_FN_
    {
        if constexpr (device_atomic)
        {
            return atomicCAS(addr, compare, value);
        }
        else if constexpr (atomic)
        {
            // like atomicCAS: returns the value at addr before the call, == compare on success
            std::atomic_ref<T>{*addr}.compare_exchange_strong(compare, value);
            return compare;
        }
        else
            static_assert(dependent_false_v<T>, "compare_and_swap requires atomic");
    }

    T& t;
};
} // namespace PHARE::core

#endif /* PHARE_CORE_OPERATORS_HPP */
