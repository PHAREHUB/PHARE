#ifndef PHARE_CORE_UTILITIES_ALLOCATORS_HPP
#define PHARE_CORE_UTILITIES_ALLOCATORS_HPP


#include <new>
#include <cstddef>
#include <cstdlib>
#include <type_traits>

#include <limits>
#include <utility>
#include <stdlib.h>   // posix_memalign
#include <sys/mman.h> // madvise

namespace PHARE::core
{

template<typename T>
class Allocator
{
    using This = Allocator<T>;

public:
    using pointer         = T*;
    using reference       = T&;
    using value_type      = T;
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;

    template<typename U>
    struct rebind
    {
        using other = Allocator<U>;
    };

    T* allocate(std::size_t const n) const
    {
        return static_cast<T*>(::operator new(n * sizeof(T)));
    }

    void deallocate(T* const p) noexcept
    {
        if (p)
            ::operator delete(p);
    }
    void deallocate(T* const p, std::size_t /*n*/) noexcept // needed from std::
    {
        deallocate(p);
    }

    bool operator!=(This const& that) const { return !(*this == that); }
    bool operator==(This const& /*that*/) const
    {
        return true; // stateless
    }
};

template<typename T>
class NonConstructingAllocator : public Allocator<T>
{
    using This = NonConstructingAllocator<T>;

public:
    template<typename U>
    struct rebind
    {
        using other = NonConstructingAllocator<U>;
    };

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

template<typename T, std::size_t huge_page_size = 4096>
class HugePageAllocator : public Allocator<T>
{
    using This = HugePageAllocator<T, huge_page_size>;

public:
    template<typename U>
    struct rebind
    {
        using other = HugePageAllocator<U>;
    };

    HugePageAllocator() = default;
    template<class U>
    constexpr HugePageAllocator(HugePageAllocator<U> const&) noexcept
    {
    }

    T* allocate(std::size_t n)
    {
#ifndef MADV_HUGEPAGE
        throw std::runtime_error("HugePageAllocator not available");
#endif

        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::bad_alloc();
        void* p = nullptr;

        if (posix_memalign(&p, huge_page_size, n * sizeof(T)) != 0)
            throw std::bad_alloc{};

#ifdef MADV_HUGEPAGE
        madvise(p, n * sizeof(T), MADV_HUGEPAGE);
#endif // MADV_HUGEPAGE
        if (p == nullptr)
            throw std::bad_alloc();
        return static_cast<T*>(p);
    }

    void deallocate(T* const p) noexcept
    {
        if (p)
            std::free(p);
    }
    void deallocate(T* const p, std::size_t /*n*/) noexcept // needed from std::
    {
        deallocate(p);
    }

    bool operator!=(This const& that) const { return !(*this == that); }
    bool operator==(This const& /*that*/) const
    {
        return true; // stateless
    }
};

template<typename T, std::size_t huge_page_size = 4096>
class NonConstructingHugePageAllocator : public NonConstructingAllocator<T>
{
    using This = NonConstructingHugePageAllocator<T, huge_page_size>;

public:
    template<typename U>
    struct rebind
    {
        using other = NonConstructingHugePageAllocator<U>;
    };

    NonConstructingHugePageAllocator() = default;
    template<class U>
    constexpr NonConstructingHugePageAllocator(NonConstructingHugePageAllocator<U> const&) noexcept
    {
    }

    T* allocate(std::size_t n)
    {
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::bad_alloc();
        void* p = nullptr;
        if (posix_memalign(&p, huge_page_size, n * sizeof(T)) != 0)
            throw std::bad_alloc{};

#ifdef MADV_HUGEPAGE
        madvise(p, n * sizeof(T), MADV_HUGEPAGE);
#endif // MADV_HUGEPAGE
        if (p == nullptr)
            throw std::bad_alloc();
        return static_cast<T*>(p);
    }

    void deallocate(T* const p) noexcept
    {
        if (p)
            std::free(p);
    }
    void deallocate(T* const p, std::size_t /*n*/) noexcept // needed from std::
    {
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


} // namespace PHARE::core


#endif /*PHARE_CORE_UTILITIES_ALLOCATORS_HPP*/