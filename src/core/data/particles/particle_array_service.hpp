#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERVICE_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERVICE_HPP

#include <array>
#include <cstdint>
#include <tuple>

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_detail.hpp"
#include "core/vector.hpp"

namespace PHARE::core
{


struct ParticleArrayService
{
    // optionally implemented on arrays
    template<std::uint8_t Phase, auto type, typename ParticleArray_t>
    auto _PHARE_ALL_FN_ static sync(ParticleArray_t& p)
    {
        PHARE_LOG_SCOPE(3, "ParticleArrayService::sync");

        if constexpr (is_tiled(ParticleArray_t::layout_mode))
            p.template sync<Phase, type>();
    }

    template<auto type, typename ParticleArray_t>
    auto static reserve_ppc_in(ParticleArray_t& p, std::size_t const ppc)
    {
        PHARE_LOG_SCOPE(3, "ParticleArrayService::reserve_ppc_in");

        using enum LayoutMode;
        if constexpr (any_in(ParticleArray_t::layout_mode, AoSPC, SoAPC)
                      or is_tiled(ParticleArray_t::layout_mode))
            p.template reserve_ppc<type>(ppc);
    }




    // noops below


    // template<std::uint8_t Phase, auto type, typename... Args>
    // void _PHARE_ALL_FN_ static sync(Args&&...)
    // {
    //     PHARE_LOG_SCOPE(3, "ParticleArrayService::sync noop");
    // }


    // template<auto type, typename... Args>
    // void static reserve_ppc_in(Args&&... args)
    // {
    //     PHARE_LOG_SCOPE(3, "ParticleArrayService::reserve_ppc_in noop");

    //     using Tuple_t         = std::tuple<Args...>;
    //     using ParticleArray_t = std::decay_t<std::tuple_element_t<0, Tuple_t>>;

    //     static_assert(not any_in(ParticleArray_t::layout_mode, AoSTS));
    // }


    //
};



} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERVICE_HPP*/
