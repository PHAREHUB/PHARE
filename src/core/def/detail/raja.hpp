#ifndef PHARE_CORE_DEF_DETAIL_RAJA_HPP
#define PHARE_CORE_DEF_DETAIL_RAJA_HPP

#if !defined(PHARE_HAVE_RAJA)
#if __has_include("RAJA/RAJA.hpp")
#define PHARE_HAVE_RAJA 1
#endif // __has_include("RAJA/RAJA.hpp")
#endif // !defined(PHARE_HAVE_RAJA)


#if !defined(PHARE_HAVE_RAJA)
#define PHARE_HAVE_RAJA 0
#endif

#if PHARE_HAVE_RAJA

#if defined(HAVE_SYS_TIMES_H)
#undef HAVE_SYS_TIMES_H // https://github.com/LLNL/SAMRAI/issues/93
#undef HAVE_UNISTD_H
#endif

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Collectives.h" // tbox::parallel_synchronize();
#include "RAJA/RAJA.hpp"
#include "RAJA/config.hpp"

#define PHARE_WITH_RAJA(...) __VA_ARGS__
#else

#define PHARE_WITH_RAJA(...)
#endif


#if PHARE_HAVE_RAJA

namespace PHARE::core::raja
{
inline void sync()
{
    SAMRAI::tbox::parallel_synchronize();
}

#if !defined(RAJA_ENABLE_CUDA)

template<typename T>
void set(T* const dst, T val, std::size_t const size)
{
    std::fill(dst, dst + size, val);
}

template<typename T>
void copy(T* const dst, T const* const src, std::size_t const size)
{
    std::copy(src, src + size, dst);
}

#else  // defined(RAJA_ENABLE_CUDA)
template<typename T>
void set(T* const dst, T val, std::size_t size)
{
    std::vector<T> src(size, val); // :(
    // RAJA::resources::Cuda{}.memset(dst, val, sizeof(T) * size); // is int only
    RAJA::resources::Cuda{}.memcpy(dst, src.data(), sizeof(T) * size);
    sync();
}

template<std::size_t N = 1024, typename Fn>
void exec(Fn& fn, RAJA::resources::Cuda& res, std::size_t size)
{
    using EXEC_POLICY = RAJA::cuda_exec_async<N>;
    RAJA::forall<EXEC_POLICY>(res, RAJA::RangeSegment(0, size),
                              [=] RAJA_DEVICE(int i) mutable { fn(i); });
    sync();
}

template<std::size_t N = 1024, typename Fn>
void exec(Fn&& fn, RAJA::resources::Cuda& res, std::size_t const size)
{
    exec<N>(fn, res, size);
}

template<std::size_t N = 1024, typename Fn>
void exec(Fn&& fn, std::size_t size)
{
    RAJA::resources::Cuda res;
    exec<N>(fn, res, size);
}


template<typename T>
void copy(T* const dst, T const* const src, std::size_t const size)
{
    RAJA::resources::Cuda{}.memcpy(dst, src, sizeof(T) * size); // is async
    sync();
}


template<typename CudaRes, typename Type>
auto allocate_copy(CudaRes& res, Type* src, std::size_t size = 1)
{
    assert(size);
    Type* dst = res.template allocate<Type>(size);
    res.memcpy(dst, src, sizeof(Type) * size);
    return dst;
}


template<typename CudaRes, typename Type>
auto allocate_copy(CudaRes& res, std::vector<Type>& vec)
{
    return allocate_copy(res, vec.data(), vec.size());
}

template<typename CudaRes, typename... Types>
auto allocate_copy(CudaRes& res, std::tuple<Types...>& tup)
{
    return allocate_copy(res, &tup);
}

template<typename CudaRes, typename... Ptrs>
auto deallocate(CudaRes& res, Ptrs&... ptrs)
{
    (res.deallocate(ptrs), ...);
}
#endif // defined(RAJA_ENABLE_CUDA)

} // namespace PHARE::core::raja


#endif /* PHARE_HAVE_RAJA */
#endif /* PHARE_CORE_DEF_DETAIL_RAJA_HPP */
