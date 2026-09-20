#ifndef PHARE_CORE_OPERATORS_HPP
#define PHARE_CORE_OPERATORS_HPP

#include "core/utilities/types.hpp"

#include <atomic>

namespace PHARE::core
{
template<typename T, bool atomic>
struct Operators
{
    T static constexpr ONE = 1;

    static_assert(not std::is_const_v<T>); // doesn't make sense

    void operator+=(T const& v)
    {
        if constexpr (atomic)
        {
            auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(&t);
            T tmp          = atomic_t.load();
            while (!atomic_t.compare_exchange_weak(tmp, tmp + v))
            {
            }
        }
        else
            t += v;
    }
    void operator+=(T const&& v) { (*this) += v; }

    void operator-=(T const& v)
    {
        if constexpr (atomic)
        {
            auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(&t);
            T tmp          = atomic_t.load();
            while (!atomic_t.compare_exchange_weak(tmp, tmp - v))
            {
            }
        }
        else
            t -= v;
    }
    void operator-=(T const&& v) { (*this) -= v; }

    auto increment_return_old() // postfix increment
    {
        if constexpr (atomic)
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

    auto static compare_and_swap(T* addr, T compare, T value)
    {
        if constexpr (atomic)
        {
            // returns the value at addr before the call, == compare on success
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
