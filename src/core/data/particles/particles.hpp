#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLES_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLES_HPP


// #include "core/data/particles/particle.hpp"
// #include "core/data/particles/particle_array.hpp"
// #include "core/data/particles/particle_array_soa.hpp"
// #include "core/data/particles/particle_packer.hpp"

// namespace PHARE::core
// {
// template<std::size_t dim, typename T>
// inline constexpr auto is_phare_particle_type = //
//     std::is_same_v<Particle<dim>, T> or        //
//     is_phare_particle_view_type<dim, T>;


// template<std::size_t dim, template<std::size_t> typename ParticleA,
//          template<std::size_t> typename ParticleB>
// typename std::enable_if_t<
//     is_phare_particle_type<dim, ParticleA<dim>> and is_phare_particle_type<dim, ParticleB<dim>>,
//     bool>
// operator==(ParticleA<dim> const& particleA, ParticleB<dim> const& particleB)
// {
//     return particleA.weight == particleB.weight and //
//            particleA.charge == particleB.charge and //
//            particleA.iCell == particleB.iCell and   //
//            particleA.delta == particleB.delta and   //
//            particleA.v == particleB.v;
// }

// } // namespace PHARE::core


// namespace std
// {
// using namespace PHARE::core;

// template<std::size_t dim, bool _const_>
// Particle<dim> copy(ParticleViewBase<dim, _const_> const& from)
// {
//     return {from.weight, from.charge, from.iCell, from.delta, from.v};
// }

// template<std::size_t dim, template<std::size_t> typename Particle_t>
// Particle<dim> copy(Particle<dim> const& from)
// {
//     return {from.weight, from.charge, from.iCell, from.delta, from.v};
// }

// } // namespace std

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLES_HPP */
