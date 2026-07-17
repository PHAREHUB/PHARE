#ifndef PHARE_CORE_OPERATORS_HPP
#define PHARE_CORE_OPERATORS_HPP

#include "core/def.hpp"
#include "core/utilities/span.hpp"

#ifndef PHARE_HAVE_GPU
#define PHARE_HAVE_GPU 0
#endif

#include <atomic>


namespace PHARE::core
{
template<typename T, bool atomic, bool GPU>
struct Operators
{
    T static constexpr ONE = 1;

    static_assert(not std::is_const_v<T>); // doesn't make sense

    void operator+=(T const& v) _PHARE_ALL_FN_
    {
        if constexpr (GPU and atomic)
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
        if constexpr (GPU and atomic)
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
    void operator-=(T const&& v) _PHARE_ALL_FN_ { (*this) += v; }

    auto increment_return_old() _PHARE_ALL_FN_ // postfix increment
    {
        if constexpr (GPU and atomic)
        {
            auto o = atomicAdd(&t, ONE);
            PHARE_ASSERT(o < t);
            return o;
        }
        else if constexpr (atomic)
        {
            // update to C++20 atomic_ref (same as below but less hacky)
            auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(&t); // UB || !UB ?
            T tmp          = atomic_t.load();
            while (!atomic_t.compare_exchange_weak(tmp, tmp + 1)) {}
            return tmp;
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
        if constexpr (GPU and atomic)
        {
            return atomicCAS(addr, compare, value);
        }
        else if constexpr (atomic)
        {
            auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(addr);
            T old          = atomic_t.load();
            if (atomic_t.compare_exchange_weak(compare, value))
                return old;
            return value;
        }
    }

    T& t;
};
} // namespace PHARE::core

#endif /* PHARE_CORE_OPERATORS_HPP */
