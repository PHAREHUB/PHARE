#ifndef PHARE_CORE_UTILITIES_PARTITIONNER_PARTITIONNER_H
#define PHARE_CORE_UTILITIES_PARTITIONNER_PARTITIONNER_H

#include <vector>

#include "core/utilities/meta/meta_utilities.h"
#include "core/utilities/particle_selector/particle_selector.h"

namespace PHARE
{
namespace core
{
    /** A Partitionner taks a range of Particles and a container of boxes
     * and partitions the range so that particles are sorted according to the
     * box they belong to. Particles not belonging to any box will be
     * between the last iterator and end.
     *
     * Example with {Box1, Box2, Box3}, assuming 25, 34 and 72 particles are in each
     * box, respectively, and 21 particles are in none. In this case the function returns:
     *
     * {begin, pivot1, pivot2, pivot3}
     *
     * particles in Box1 are in [begin, pivot1[
     * particles in Box2 are in [pivot1, pivot2[
     * particles in Box3 are in [pivot2, pivot3[
     * particles in none of the boxes are in [pivot3, end[
     *
     * this function is useful after pushing the particles, when the given range
     * [begin, end[ are leaving particles, some are in physical boundary boxes, some
     * are leaving the patch but not through physical boundaries.
     *
     */
    template<typename ParticleIterator, typename BoxContainer,
             is_iterable<BoxContainer> = dummy::value>
    auto partitionner(ParticleIterator begin, ParticleIterator end, BoxContainer boxes)
    {
        std::vector<ParticleIterator> iterators;
        iterators.push_back(begin);
        auto pivot = begin;

        for (auto const& box : boxes)
        {
            pivot = std::partition(pivot, end, makeSelector(box));
            iterators.push_back(pivot);
        }

        return iterators;
    }
} // namespace core
} // namespace PHARE
#endif
