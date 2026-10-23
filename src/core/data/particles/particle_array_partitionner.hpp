#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PARTITIONNER_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PARTITIONNER_HPP

#include "core/utilities/box/box.hpp"
#include "core/utilities/partitionner/partitionner.hpp"

namespace PHARE::core
{
/** eg
  auto iterators = partition(particles, box);
  iterators[0] = in domain
  iterators[1] = ghost layer
  iterators[1].end() to particles.end() = outside ghost layer  */
template<std::size_t particle_ghost_width, typename ParticleArray, typename Box_t>
auto static partition(ParticleArray& particles, Box_t const& box)
{
    std::vector<Box_t> partition_boxes{/*domain*/ box,
                                       /*ghostLayer*/ grow(box, particle_ghost_width)};
    return partitionner(particles.begin(), particles.end(), partition_boxes);
}

/** eg
  auto iterators = partition(particles, box);
  iterators[0] = inner domain - no neighbour ghostbox overlap
  iterators[...] = contigous memory for copying to other patches - not sure how to use yet
  iterators.back().end() to particles.end() = outside domain */
template<std::size_t particle_ghost_width, typename ParticleArray, typename Box_t>
auto static partition(ParticleArray& particles, Box_t const& box, std::vector<Box_t> neighbor_boxes)
{
    auto neighbor_ghost_boxes
        = generate_from([](auto box) { return box.grow(particle_ghost_width); }, neighbor_boxes);
    PHARE_ASSERT(all_overlaps(neighbor_ghost_boxes, box));
    auto overlaps = distinct_overlaps(neighbor_ghost_boxes, box);

    std::vector<Box_t> partition_boxes{shrink(box, particle_ghost_width)};
    std::copy(overlaps.begin(), overlaps.end(), std::back_inserter(partition_boxes));

    return partitionner(particles.begin(), particles.end(), partition_boxes);
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PARTITIONNER_HPP */
