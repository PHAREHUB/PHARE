#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP

#include <cstddef>

#include "core/data/particles/particle_array_aos.hpp"
#include "core/data/particles/particle_array_soa.hpp"

namespace PHARE::core
{
template<std::size_t dim, bool SOA>
using ParticleArrayBase
    = std::conditional_t<SOA, SoAParticles<SoAVector<dim>>, AoSParticles<AoSVector<dim>>>;

template<std::size_t dim, bool SOA = false>
class ParticleArray : public ParticleArrayBase<dim, SOA>
{
public:
    using This  = ParticleArray<dim, SOA>;
    using Super = ParticleArrayBase<dim, SOA>;

    static constexpr bool is_contiguous = Super::is_contiguous;
    static constexpr auto dimension     = Super::dimension;

    using typename Super::value_type;

    template<std::size_t size>
    using array_type = typename Super::template array_type<size>;

    using typename Super::const_iterator;
    using typename Super::iterator;

    // using Super::begin;
    // using Super::end;

    ParticleArray() {}

    ParticleArray(std::size_t size)
        : Super{size}
    {
    }

    template<typename Particle_t>
    ParticleArray(std::size_t size, Particle_t&& particle)
        : Super{size, std::forward<Particle_t>(particle)}
    {
    }

    template<typename It>
    ParticleArray(It start, It end)
        : Super{start, end}
    {
    }
};


} // namespace PHARE::core


namespace PHARE::core
{
template<std::size_t dim, bool SOA>
void empty(ParticleArray<dim, SOA>& array)
{
    array.clear();
}

template<std::size_t dim, bool SOA>
void swap(ParticleArray<dim, SOA>& array1, ParticleArray<dim, SOA>& array2)
{
    array1.swap(array2);
}

} // namespace PHARE::core




#endif
