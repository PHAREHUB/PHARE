#ifndef PHARE_CORE_UTILITIES_MEMORY_HPP
#define PHARE_CORE_UTILITIES_MEMORY_HPP


#include "core/def/kokkos_tools.hpp"
#include "core/def/phare_config.hpp"



namespace PHARE::core::gpu
{


template<typename T0, typename T1, typename Size>
void copy(T0* dst, T1* src, Size const size)
{
    if constexpr (CompileOptions::WithRAJA)
    {
        PHARE_WITH_RAJA(PHARE::core::raja::copy(dst, src, size));
    }
    else if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU({
            assert(
                not(mkn::gpu::Pointer{dst}.is_host_ptr() and mkn::gpu::Pointer{src}.is_host_ptr()));
            mkn::gpu::copy(dst, src, size);
        })
    }
    else
        throw std::runtime_error("Vector::copy NO ALTERNATIVE");
}


template<typename Vector, typename Value>
void fill(Vector& dst, Value const val)
{
    if constexpr (CompileOptions::WithRAJA)
    {
        PHARE_WITH_RAJA(PHARE::core::raja::set(dst.data(), val, dst.size()));
    }
    else if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU(mkn::gpu::fill(dst, val));
    }
    else if (CompileOptions::WithKokkos)
    {
        PHARE_WITH_KOKKOS({ //
            auto data = dst.data();
            Kokkos::parallel_for("fill", 5, KOKKOS_LAMBDA(int i) { data[i] = val; });
        })
    }
    else
        throw std::runtime_error("Vector::fill NO ALTERNATIVE");
}



} // namespace PHARE::core::gpu


namespace PHARE::core::mem
{

template<auto alloc_mode, typename T0, typename T1, typename Size>
void copy(T0* dst, T1* src, Size const size)
{
    if constexpr (any_in(AllocatorMode::GPU_UNIFIED, alloc_mode))
        PHARE::core::gpu::copy(dst, src, size);
    else
        std::copy(src, src + size, dst);
}

template<auto alloc_mode0, auto alloc_mode1, typename T0, typename T1, typename Size>
void copy(T0* dst, T1* src, Size const size)
{
    if constexpr (any_in(AllocatorMode::GPU_UNIFIED, alloc_mode0, alloc_mode1))
        PHARE::core::gpu::copy(dst, src, size);
    else
        std::copy(src, src + size, dst);
}

} // namespace PHARE::core::mem


#if PHARE_HAVE_KOKKOS

namespace PHARE
{

template<typename T>
class KokkosAllocator
{
    using This = KokkosAllocator<T>;

public:
    using pointer         = T*;
    using reference       = T&;
    using value_type      = T;
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;

    template<typename U>
    struct rebind
    {
        using other = KokkosAllocator<U>;
    };

    T* allocate(std::size_t const n) const
    {
        if (n == 0)
            return nullptr;

        void* ptr = Kokkos::kokkos_malloc(sizeof(T) * n);

        if (!ptr)
            throw std::bad_alloc();
        return static_cast<T*>(ptr);
    }

    void deallocate(T* const p) noexcept
    {
        if (p)
            Kokkos::kokkos_free(p);
    }
    void deallocate(T* const p, std::size_t /*n*/) noexcept
    { // needed from std::
        deallocate(p);
    }

    template<typename U, typename... Args>
    void construct(U* ptr, Args&&... args)
    {
        ::new ((void*)ptr) U(std::forward<Args>(args)...);
    }
    template<typename U>
    void construct(U* /*ptr*/) noexcept(std::is_nothrow_default_constructible<U>::value)
    {
    }

    bool operator!=(This const& that) const { return !(*this == that); }

    bool operator==(This const& /*that*/) const
    {
        return true; // stateless
    }
};




} // namespace PHARE


#endif


#endif /* PHARE_CORE_UTILITIES_MEMORY_HPP */
