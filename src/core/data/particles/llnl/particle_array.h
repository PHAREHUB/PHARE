#ifndef PHARE_CORE_LLNL_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_LLNL_PARTICLES_PARTICLE_ARRAY_H

#ifndef HAVE_UMPIRE
#error // expected
#endif

#include <atomic>
#include <cstddef>

#include "core/def.hpp"
#include "core/data/particles/particle_array.hpp"

namespace PHARE::core::llnl
{
template<typename Particle>
struct ABufferedParticleVector
{
    using Allocator = umpire::TypedAllocator<Particle>;

    auto get_allocator()
    {
        auto& rm = umpire::ResourceManager::getInstance();
        assert(rm.isAllocator("PHARE::data_allocator"));
        return Allocator{rm.getAllocator("PHARE::data_allocator")};
    }

    ABufferedParticleVector(std::size_t size_, double buffer_by_, double realloc_by_)
        : size{size_}
        , n_elements{size_}
        , buffer_by{realloc_by_}
        , realloc_by{realloc_by_}
        , capacity{static_cast<std::size_t>(size + size * buffer_by)}
        , particles{get_allocator()}
    {
        if (size_)
            particles.resize(capacity);

        assert(particles.size() == capacity);
        assert(buffer_by < 1 and buffer_by > 0);
        assert(realloc_by < 1 and realloc_by > 0);
    }


    ABufferedParticleVector(double buffer_by, double realloc_by)
        : ABufferedParticleVector{0, buffer_by, realloc_by}
    {
    }

    bool check()
    {
        assert(this);
        bool realloc_more = n_elements >= size + size * realloc_by;
        KLOG(INF) << n_elements;
        KLOG(INF) << realloc_more;

        if (realloc_more)
        {
            size = n_elements;
            capacity += size * buffer_by;
            particles.reserve(capacity);
        }

        // bool realloc_less = false;
        // if (realloc_less)
        // {
        //     // TODO?
        // }
    }

    void clear(bool force = false)
    {
        n_elements = 0;
        if (force)
            particles.clear();
    }

    void erase(std::size_t index) __device__
    {
        assert(false);
        // auto top = atomicSub(&info[0], 1);
        // if (index != top)
        //     particles[index] = top;
    }
    void reserve(std::size_t size) { particles.reserve(size); }


    void push_back(Particle const& particle) __device__
    {
        assert(false);
        // auto index       = atomicAdd(&info[0], 1);
        // particles[index] = particle;
    }
    void push_back(Particle&& particle) __device__ { push_back(particle); }
    void push_back(Particle const& particle) __host__
    {
        // check(); // this can only work with CUDA unified memory
        // particles[++n_elements] = particle;
        KLOG(INF) << particles.size();
        particles.push_back(particle);
        KLOG(INF) << particles.size();
    }
    void push_back(Particle&& particle) __host__ { push_back(particle); }


    std::size_t size, n_elements;
    double buffer_by, realloc_by;
    std::size_t capacity;
    std::vector<Particle, Allocator> particles;
};


template<typename Particle>
struct BufferedParticleVector : public ABufferedParticleVector<Particle>
{
    using Super    = ABufferedParticleVector<Particle>;
    using iterator = Particle*;
    using Super::particles;


    BufferedParticleVector(std::size_t size, double buffer_by = .2, double realloc_by = .1) __host__
        : ABufferedParticleVector<Particle>{size, buffer_by, realloc_by}
    {
    }

    BufferedParticleVector(double buffer_by = .2, double realloc_by = .1) __host__
        : ABufferedParticleVector<Particle>{buffer_by, realloc_by}
    {
    }

    auto& size() const __device__ __host__ { return Super::size; }
    auto& operator[](std::size_t i) const __device__ { return particles[i]; }
    auto& operator[](std::size_t i) __device__ { return particles[i]; }
};

template<typename Particle, typename Vector_ = BufferedParticleVector<Particle>>
class ParticleArray
{
    bool def_const = 0; // TORM

public:
    static constexpr bool is_host_mem   = false;
    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = Particle::dimension;
    using Particle_t                    = Particle;
    using Vector                        = Vector_;
    using iterator                      = typename Vector::iterator;
    using value_type                    = Particle_t;


    ParticleArray()
        : vector{std::make_unique<Vector>()}
    {
        KLOG(INF);
        def_const = 1;
        this->vector->reserve(100);
    }

    ParticleArray(std::size_t size)
        : vector{std::make_unique<Vector>(size)}
    {
        KLOG(INF) << size;
    }

    ParticleArray(ParticleArray const& that)               = delete;
    ParticleArray(ParticleArray&& that)                    = delete;
    ParticleArray(std::size_t size, Particle_t&& particle) = delete;


    std::size_t size() const
    {
        KLOG(INF) << def_const;
        check();
        return vector->size();
    }

    void clear() { vector->clear(); }

    auto& operator=(std::vector<Particle_t>&& input)
    {
        // this->particles = std::move(vector);
        this->vector = std::make_unique<Vector>(input.size());
        KLOG(INF) << input.size();
        KLOG(INF) << vector->size();

        KLOG(INF) << &input;
        KLOG(INF) << &vector->particles;
        KLOG(INF) << vector->particles.data();

        assert(vector != nullptr);
        assert(input.size() == vector->size());

        RAJA::resources::Cuda{}.memcpy(
            /*device pointer*/ vector->particles.data(),
            /*host pointer*/ input.data(),
            /*size in bytes*/ sizeof(Particle_t) * input.size());

        return *this;
    }

    operator bool() const { return vector != nullptr; }
    void check() const { assert(bool{*this}); }

    // auto& operator[](std::size_t i) _PHARE_ALL_FN_
    // {
    //     check();
    //     return vector[i];
    // }
    // auto& operator[](std::size_t i) const _PHARE_ALL_FN_
    // {
    //     check();
    //     return vector[i];
    // }

    auto data() const
    {
        KLOG(INF);
        check();
        return vector->particles.data();
    }
    auto data()
    {
        KLOG(INF);
        check();
        return vector->particles.data();
    }
    void erase(std::size_t index)
    {
        KLOG(INF);
        check();
        vector->erase(index);
    }

    void push_back(Particle&& p)
    {
        KLOG(INF);
        check();
        vector->push_back(p);
    }
    void push_back(Particle const& p)
    {
        KLOG(INF);
        check();
        vector->push_back(p);
    }

    void swap(ParticleArray<Particle>& that)
    {
        KLOG(INF);
        check();
        this->vector->particles.swap(that.vector->particles);
    }

    auto operator()() const
    {
        KLOG(INF) << size();
        check();
        core::ParticleArray<Particle> particles(size());
        PHARE_WITH_RAJA(PHARE::raja::copy(particles.data(), vector->particles.data(), size()));
        return particles;
    }



    auto begin() const
    {
        check();
        return vector->particles.begin();
    }
    auto begin()
    {
        check();
        return vector->particles.begin();
    }

    auto end() const
    {
        check();
        return vector->particles.end();
    }
    auto end()
    {
        check();
        return vector->particles.end();
    }

    auto& operator[](std::size_t i) const
    {
        check();
        return vector->particles[i];
    }
    auto& operator[](std::size_t i)
    {
        check();
        return vector->particles[i];
    }

private:
    std::unique_ptr<Vector> vector;
};

} // namespace PHARE::core::llnl


namespace PHARE::core
{
template<typename Particle>
void empty(llnl::ParticleArray<Particle>& array)
{
    // assert(false);
    array.clear();
}

template<typename Particle>
void swap(llnl::ParticleArray<Particle>& array1, llnl::ParticleArray<Particle>& array2)
{
    // assert(false);
    array1.swap(array2);
}

template<typename Particle>
void append(llnl::ParticleArray<Particle> const& src, llnl::ParticleArray<Particle>& dst)
{
    std::copy(std::begin(src), std::end(src), std::back_inserter(dst));
}

} // namespace PHARE::core

#endif /*PHARE_CORE_LLNL_PARTICLES_PARTICLE_ARRAY_H*/
