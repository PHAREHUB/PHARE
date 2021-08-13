#ifndef PHARE_CORE_OPERATORS_H
#define PHARE_CORE_OPERATORS_H

#include "core/def.h"

#if defined(PHARE_WITH_GPU)

#include <cuda_runtime.h> // TORM

namespace PHARE
{
struct GPU_Operators
{
    template<typename T>
    struct Operator
    {
        void operator+=(T const& v) _PHARE_FN_SIG_ { t += v; }

        T& t;
    };

    template<typename T>
    struct AtomicOperator
    {
        AtomicOperator() _PHARE_GPU_FN_DEV_ = default;
        AtomicOperator() _PHARE_GPU_FN_HST_ = delete;

        void operator+=(T const& v) _PHARE_GPU_FN_DEV_ { atomicAdd(&t, v); }
        void operator+=(T const&& v) _PHARE_GPU_FN_DEV_ { atomicAdd(&t, v); }

        void operator+=(T const& v) _PHARE_GPU_FN_HST_
        {
            throw_runtime_error("1 PHARE NOT BUILT WITH PHARE_WITH_GPU");
        }
        void operator+=(T const&& v) _PHARE_GPU_FN_HST_
        {
            throw_runtime_error("2 PHARE NOT BUILT WITH PHARE_WITH_GPU");
        }

        T& t;
    };
};


} // namespace PHARE

#else  // ! PHARE_WITH_GPU

namespace PHARE
{
struct GPU_Operators
{
    template<typename T>
    struct Operator
    {
        void operator+=(T const& v)
        {
            throw_runtime_error("3 PHARE NOT BUILT WITH PHARE_WITH_GPU");
        }
    };

    template<typename T>
    struct AtomicOperator
    {
        void operator+=(T const& v)
        {
            throw_runtime_error("4 PHARE NOT BUILT WITH PHARE_WITH_GPU");
        }
        void operator+=(T const&& v)
        {
            throw_runtime_error("5 PHARE NOT BUILT WITH PHARE_WITH_GPU");
        }
    };
};
} // namespace PHARE
#endif // PHARE_WITH_GPU

struct CPU_Operators
{
    template<typename T>
    struct Operator
    {
        void operator+=(T const& v) { t += v; }
        void operator+=(T const&& v) { t += v; }
        // void operator-(T& v) { t -= v; }

        T& t;
    };

    template<typename T>
    struct AtomicOperator
    {
        void operator+=(T& v) { __sync_fetch_and_add(&t, v); }

        T& t;
    };
};

#endif /*PHARE_CORE_OPERATORS_H*/
