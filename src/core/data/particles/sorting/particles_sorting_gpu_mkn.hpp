#ifndef PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_GPU_HPP
#define PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_GPU_HPP

#if PHARE_HAVE_THRUST

#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <thrust/sort.h>
#include <thrust/gather.h>
#include <thrust/sequence.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>


namespace PHARE::core::detail
{


template<typename ParticleArray>
class ParticleSorter<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray>
{
public:
    using box_t = Box<int, ParticleArray::dimension>;

    void operator()(std::uint64_t l, std::uint64_t r)
    {
        PHARE_LOG_SCOPE(1, "ParticleSorter<AllocatorMode::GPU_UNIFIED>");

        if constexpr (ParticleArray::layout_mode == LayoutMode::AoS)
        {
            auto ps = particles.view();
            thrust::sort(thrust::device, ps.begin() + l, ps.begin() + r,
                         [cf = cell_flattener] _PHARE_ALL_FN_(auto const& a, auto const& b) {
                             return cf(a.iCell()) < cf(b.iCell());
                         });
        }
        else if constexpr (ParticleArray::layout_mode == LayoutMode::SoA)
        {
            auto it = SoAIteratorAdaptor::make(particles);
            thrust::sort(thrust::device, it + l, it + r,
                         [cf = cell_flattener] _PHARE_ALL_FN_(auto const& a, auto const& b) {
                             return cf(SoAIteratorAdaptor::iCell(a))
                                    < cf(SoAIteratorAdaptor::iCell(b));
                         });
        }
        else
            throw std::runtime_error("no impl");
    }

    void operator()() { (*this)(0, particles.size()); }


    auto static constexpr by_deltas() _PHARE_ALL_FN_
    {
        return [](auto const& a, auto const& b) -> bool {
            return as_tuple(a.delta()) < as_tuple(b.delta());
        };
    }
    void by_deltas(std::uint64_t l, std::uint64_t r)
    {
        PHARE_LOG_SCOPE(1, "ParticleSorter<AllocatorMode::GPU_UNIFIED>");
        if constexpr (ParticleArray::layout_mode == LayoutMode::AoS)
        {
            auto ps = particles.view();
            thrust::sort(thrust::device, ps.begin() + l, ps.begin() + r, by_deltas());
        }
        else if constexpr (ParticleArray::layout_mode == LayoutMode::SoA)
        {
            auto it = SoAIteratorAdaptor::make(particles);
            thrust::sort(thrust::device, it + l, it + r, [](auto const& a, auto const& b) -> bool {
                return as_tuple(SoAIteratorAdaptor::delta(a))
                       < as_tuple(SoAIteratorAdaptor::delta(b));
            });
        }
        else
            throw std::runtime_error("no impl");
    }

    ParticleSorter& by_delta()
    {
        PHARE_LOG_SCOPE(1, "ParticleSorter<AllocatorMode::GPU_UNIFIED>");
        // assumes already sorted by icell
        if (particles.size() == 0)
            return *this;

        using enum LayoutMode;

        if constexpr (any_in(ParticleArray::layout_mode, AoSPC))
            for (auto const& bix : particles.local_box())
            {
                auto& ps = particles(bix);
                thrust::sort(thrust::device, ps.data(), ps.data() + ps.size(), by_deltas());
            }
        else if constexpr (any_in(ParticleArray::layout_mode, SoAPC))
            for (auto const& bix : particles.local_box())
            {
                auto& ps = particles(bix);
                auto it  = SoAIteratorAdaptor::make(ps);
                thrust::sort(thrust::device, it, it + ps.size(),
                             [](auto const& a, auto const& b) -> bool {
                                 return as_tuple(SoAIteratorAdaptor::delta(a))
                                        < as_tuple(SoAIteratorAdaptor::delta(b));
                             });
            }
        else if constexpr (any_in(ParticleArray::layout_mode, AoSTS, SoATS))
        {
            throw std::runtime_error("no impl");
        }
        else
        {
            auto const end = particles.end();
            auto beg       = particles.begin();
            auto lst       = particles.begin();
            while (lst != end)
            {
                lst = beg + 1;
                while (lst != end and lst.iCell() == beg.iCell())
                    ++lst;
                by_deltas(it_dist(particles.begin(), beg), it_dist(particles.begin(), lst));
                beg = lst;
            }
        }

        return *this;
    }


    void sort_by_key(std::uint64_t l, std::uint64_t r) // unused but for reference
    {
        auto ps = particles.view();
        std::vector<int, mkn::gpu::ManagedAllocator<int>> flats(r - l);
        auto fv = flats.data();
        mkn::gpu::GDLauncher{flats.size()}([=, cf = cell_flattener] _PHARE_ALL_FN_() {
            auto idx            = mkn::gpu::idx() + l;
            fv[mkn::gpu::idx()] = cf(ps.iCell(idx));
        });
        thrust::device_vector<int> indices(flats.size());
        thrust::sequence(indices.begin(), indices.end());
        thrust::sort_by_key(thrust::device, fv, fv + flats.size(), indices.begin());
        thrust::gather(thrust::device, indices.begin(), indices.end(), ps.begin() + l,
                       ps.begin() + l);
    }


    ParticleArray& particles;
    box_t& domain_box;
    CellFlattener<box_t> cell_flattener{domain_box};
};


// template<typename ParticleArray>
// class ParticleSorter<AllocatorMode::GPU, /*impl = */ 0, ParticleArray>
// {
// public:
//     void operator()(std::int64_t /*l*/, std::int64_t /*r*/) { std::abort(); }

//     void operator()() { (*this)(0, particles.size() - 1); }

//     ParticleArray& particles;
// };

} // namespace PHARE::core::detail


#endif // PHARE_HAVE_THRUST

#endif /*PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_GPU_HPP*/
