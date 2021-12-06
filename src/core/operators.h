#ifndef PHARE_CORE_OPERATORS_H
#define PHARE_CORE_OPERATORS_H

#include "core/def.h"

#include <atomic>

namespace PHARE::core
{
template<typename T, bool atomic = false>
struct Operators
{
    void operator+=(T const& v)
    {
        if constexpr (atomic)
        {
            auto& atomic_t = *reinterpret_cast<std::atomic<T>*>(&t);
            T tmp          = atomic_t.load();
            while (!atomic_t.compare_exchange_weak(tmp, tmp + v)) {}
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
            while (!atomic_t.compare_exchange_weak(tmp, tmp - v)) {}
        }
        else
            t -= v;
    }
    void operator-=(T const&& v) { (*this) += v; }

    T& t;
};
} // namespace PHARE::core

#endif /*PHARE_CORE_OPERATORS_H*/
