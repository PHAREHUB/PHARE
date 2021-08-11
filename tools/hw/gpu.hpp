#ifndef PHARE_GPU_HPP
#define PHARE_GPU_HPP

#include <cstdlib>
#include <algorithm>

#include "kul/gpu.hpp"
#include "kul/span.hpp"
#include "kul/gpu/asio.hpp"
#include "kul/gpu/tuple.hpp"

#include "tests/simulator/per_test.h"

namespace PHARE::gpu
{
template<typename PHARE_TYPES>
struct PatchState
{
    using Float      = double;
    using GridLayout = typename PHARE_TYPES::GridLayout_t;

    template<typename State>
    PatchState(GridLayout const& gridLayout, State& state, std::string id_)
        : id{id_}
        , layout{gridLayout}
    {
        for (auto& pop : state.ions)
            ions.emplace_back(&pop.domainParticles());

        auto vecF = [&](auto const& EBxyz) {
            for (std::uint8_t i = 0; i < 3; i++)
                electromag.emplace_back(EBxyz[i].data(), EBxyz[i].size());
        };

        vecF(state.electromag.E);
        vecF(state.electromag.B);
    }

    std::string id;
    GridLayout const layout;
    std::vector<kul::Span<Float const, std::uint32_t>> electromag;
    std::vector<core::ParticleArray</*Float,*/ PHARE_TYPES::dimension>*> ions;
};

template<std::uint8_t dim, bool GPU>
struct Particles : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = Particles<dim, true>;

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    Particles(std::uint32_t nbr)
        : particles{nbr}
    {
    }


    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    void add(core::ParticleArray<dim> const& array, std::uint32_t i)
    {
        particles.send(array.data(), array.size(), i);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto operator()()
    {
        return Super::template alloc<gpu_t>(particles);
    }


    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    auto __device__ operator[](std::size_t idx) const
    {
        return particles[idx];
    }

    container_t<core::Particle<dim>> particles;
};


template<typename Float, typename Span_ = kul::gpu::Span<Float const, std::uint32_t>>
struct FieldInterop
{
    using Span = Span_;

    FieldInterop() __device__ = default;

    auto& __device__ operator()(std::uint32_t i) { return span[i]; }
    auto& __device__ operator()(std::uint32_t i) const { return span[i]; }
    // auto operator()(std::uint32_t, std::uint32_t) { return 2; }
    // auto operator()(std::uint32_t, std::uint32_t, std::uint32_t) { return 3; }

    Span span;
};

template<typename Float, typename Span_ = kul::gpu::Span<Float const, std::uint32_t>>
struct VecFieldInterop
{
    using Span = Span_;

    VecFieldInterop() __device__ = default;

    template<typename Electromags>
    static VecFieldInterop<Float> __device__ E(Electromags const& em)
    {
        return {{Span{em.Ex, em.info[0]}}, {Span{em.Ey, em.info[1]}}, {Span{em.Ez, em.info[2]}}};
    }
    template<typename Electromags>
    static VecFieldInterop<Float> __device__ B(Electromags const& em)
    {
        return {{Span{em.Bx, em.info[3]}}, {Span{em.By, em.info[4]}}, {Span{em.Bz, em.info[5]}}};
    }

    auto& __device__ getComponent(PHARE::core::Component XYZ)
    {
        if (XYZ == PHARE::core::Component::X)
            return x;
        else if (XYZ == PHARE::core::Component::Y)
            return y;
        return z;
    }

    auto __device__ getComponents() { return std::forward_as_tuple(x, y, z); }
    auto __device__ getComponents() const { return std::forward_as_tuple(x, y, z); }

    FieldInterop<Float, Span> x, y, z;
};


template<typename Float>
struct EMInterop
{
    EMInterop() __device__ = default;

    template<typename Electromags>
    EMInterop(Electromags&& _em) __device__ : E{VecFieldInterop<Float>::E(_em)},
                                              B{VecFieldInterop<Float>::B(_em)}
    {
    }

    VecFieldInterop<Float> E, B;
};

template<typename PHARE_TYPES, bool GPU = false>
struct Electromag : kul::gpu::DeviceClass<GPU>
{
    using Float = double;
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = Electromag<PHARE_TYPES, true>;

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    Electromag(PatchState<PHARE_TYPES>& state)
        : Ex{state.electromag[0]}
        , Ey{state.electromag[1]}
        , Ez{state.electromag[2]}
        , Bx{state.electromag[3]}
        , By{state.electromag[4]}
        , Bz{state.electromag[5]}
    {
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto operator()()
    {
        return Super::template alloc<gpu_t>(Ex, Ey, Ez, Bx, By, Bz);
    }

    struct Interop
    {
        VecFieldInterop<Float, Float*> E, B;
    };

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    auto __device__ operator()() const
    {
        return Interop{{Ex, Ey, Ez}, {Bx, By, Bz}};
    }

    container_t<Float> Ex, Ey, Ez, Bx, By, Bz;
};


template<typename PHARE_TYPES, bool GPU>
struct Electromags : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = Electromags<PHARE_TYPES, true>;
    using Float = double;

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    static auto make_shared(std::vector<PatchState<PHARE_TYPES>> const& states)
    {
        std::uint32_t n_states = static_cast<std::uint32_t>(states.size());

        std::vector<std::uint32_t> emXYZ(6, 0), emInfo(n_states * 6);
        for (std::size_t j = 0; j < n_states; ++j)
            for (std::size_t i = 0; i < states[j].electromag.size(); ++i)
            {
                auto pos    = (j * 6) + i;
                emInfo[pos] = emXYZ[i];
                emXYZ[i] += states[j].electromag[i].size();
            }

        return std::make_shared<Electromags<PHARE_TYPES, GPU>>(states, emXYZ, emInfo);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    Electromags(std::vector<PatchState<PHARE_TYPES>> const& states,
                std::vector<std::uint32_t> const& v, std::vector<std::uint32_t> const& _info)
        : Ex{v[0]}
        , Ey{v[1]}
        , Ez{v[2]}
        , Bx{v[3]}
        , By{v[4]}
        , Bz{v[5]}
        , info{_info}
    {
        for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(states.size()); ++i)
        {
            auto pos = i * 6;
            auto& em = states[i].electromag;
            Ex.send(em[0], _info[pos + 0]);
            Ey.send(em[1], _info[pos + 1]);
            Ez.send(em[2], _info[pos + 2]);
            Bx.send(em[3], _info[pos + 3]);
            By.send(em[4], _info[pos + 4]);
            Bz.send(em[5], _info[pos + 5]);
        }
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto operator()()
    {
        return Super::template alloc<gpu_t>(Ex, Ey, Ez, Bx, By, Bz, info);
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    Electromags() __device__
    {
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    EMInterop<Float> __device__ electromag(std::uint16_t i)
    {
        auto pos = i * 6;
        Electromags em;
        em.Ex   = this->Ex + this->info[pos + 0];
        em.Ey   = this->Ey + this->info[pos + 1];
        em.Ez   = this->Ez + this->info[pos + 2];
        em.Bx   = this->Bx + this->info[pos + 3];
        em.By   = this->By + this->info[pos + 4];
        em.Bz   = this->Bz + this->info[pos + 5];
        em.info = this->info + pos;
        return EMInterop<Float>{em};
    }

    container_t<Float> Ex, Ey, Ez, Bx, By, Bz;
    container_t<std::uint32_t> info;
};


template<typename PHARE_TYPES, bool GPU>
struct GridLayouts : kul::gpu::DeviceClass<GPU>
{
    static constexpr std::uint8_t dim  = PHARE_TYPES::dimension;
    static constexpr std::uint8_t dim2 = dim * 2;

    using Float = double;
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = GridLayouts<PHARE_TYPES, true>;

    using GridLayoutImpl = typename PHARE_TYPES::YeeLayout_t;
    using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    static auto make_shared(std::vector<PatchState<PHARE_TYPES>> const& states)
    {
        std::uint32_t n_states = static_cast<std::uint32_t>(states.size());
        return std::make_shared<GridLayouts<PHARE_TYPES, GPU>>(states, n_states);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    GridLayouts(std::vector<PatchState<PHARE_TYPES>> const& states, std::uint32_t n_states)
        : layouts{n_states}
    {
        for (std::uint32_t i = 0; i < n_states; ++i)
            layouts.send(&states[i].layout, 1, i);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto operator()()
    {
        return Super::template alloc<gpu_t>(layouts);
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    auto __device__ operator[](std::size_t idx) const
    {
        return layouts[idx];
    }

    typename Super::template container_t<GridLayout> layouts;
};


template<typename PHARE_TYPES, bool GPU>
struct PatchStatePerParticle : kul::gpu::DeviceClass<GPU>
{
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = PatchStatePerParticle<PHARE_TYPES, true>;
    using Float = double;

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    PatchStatePerParticle(std::uint32_t n_patches, std::uint32_t n_particles)
        : particles{1}
        , electromags{1}
        , gridLayouts{1}
        , particlePatchStateIdx{n_particles}
        , info{n_patches + 1}
    {
        info.send(&n_particles);
    }

    template<bool gpu = GPU, std::enable_if_t<!gpu, bool> = 0>
    auto operator()()
    {
        return Super::template alloc<gpu_t>(particles, electromags, gridLayouts,
                                            particlePatchStateIdx, info);
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    auto __device__ n_particles() const
    {
        return info[0];
    }

    template<bool gpu = GPU, std::enable_if_t<gpu, bool> = 0>
    auto __device__ operator[](std::size_t idx) const
    {
        return particlePatchStateIdx[idx];
    }

    container_t<Particles<PHARE_TYPES::dimension, true>> particles;
    container_t<Electromags<PHARE_TYPES, true>> electromags;
    container_t<GridLayouts<PHARE_TYPES, true>> gridLayouts;

    container_t<std::uint16_t> particlePatchStateIdx;
    container_t<std::uint32_t> info;
};

template<typename PHARE_TYPES>
struct ParticlePatchState
{
    static constexpr bool GPU    = false;
    static constexpr auto dim    = PHARE_TYPES::dimension;
    using Float                  = double;
    using GridLayout             = typename GridLayouts<PHARE_TYPES, GPU>::GridLayout;
    using GridLayouts_           = GridLayouts<PHARE_TYPES, GPU>;
    using Electromags_           = Electromags<PHARE_TYPES, GPU>;
    using Particles_             = Particles<dim, GPU>;
    using PatchStatePerParticle_ = PatchStatePerParticle<PHARE_TYPES, GPU>;

    ParticlePatchState(std::vector<PatchState<PHARE_TYPES>> const& states, bool async)
    {
        std::uint32_t n_states = static_cast<std::uint32_t>(states.size());

        for (auto const& data : states)
            for (auto const& particle_array : data.ions)
                n_particles += particle_array->size();

        if (async)
            return;

        particles = std::make_shared<Particles_>(n_particles);
        pspp      = std::make_shared<PatchStatePerParticle_>(n_states, n_particles);

        for (std::uint32_t parti = 0, i = 0; i < n_states; i++)
            for (auto const* particle_array : states[i].ions)
            {
                kul::gpu::fill_n(pspp->particlePatchStateIdx + parti, particle_array->size(), i);
                particles->add(*particle_array, parti);
                parti += particle_array->size();
            }

        pspp->particles.send((*particles)());
        pspp->gridLayouts.send((*(gridLayouts = GridLayouts_::make_shared(states)))());
        pspp->electromags.send((*(electromags = Electromags_::make_shared(states)))());
    }

    auto operator()() { return (*pspp)(); }

    std::uint32_t n_particles = 0;
    std::shared_ptr<PatchStatePerParticle_> pspp;
    std::shared_ptr<Particles_> particles;
    std::shared_ptr<GridLayouts_> gridLayouts;
    std::shared_ptr<Electromags_> electromags;
};

} // namespace PHARE::gpu

#endif /*PHARE_GPU_HPP*/
