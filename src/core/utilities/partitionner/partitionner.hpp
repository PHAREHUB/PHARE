#ifndef PHARE_CORE_UTILITIES_PARTITIONNER_PARTITIONNER_HPP
#define PHARE_CORE_UTILITIES_PARTITIONNER_PARTITIONNER_HPP

#include <vector>

// #include "core/utilities/meta/meta_utilities.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
// #include "core/data/particles/particle.hpp"

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


    template<typename ParticleIterator, typename T, std::size_t dim>
    auto partitionner(ParticleIterator begin, ParticleIterator end, Box<T, dim> const& box)
    {
        return partitionner(begin, end, box,
                            [&box](auto const& part) { return isIn(cellAsPoint(part), box); });
    }


    template<typename ParticleIterator, typename T, std::size_t dim, typename Fn>
    auto partitionner(ParticleIterator begin, ParticleIterator end, Box<T, dim> const& box, Fn&& fn)
    {
        return BoxRange<Box<T, dim>, ParticleIterator>{box, begin, std::partition(begin, end, fn)};
    }

    template<typename Particles, typename T, std::size_t dim, typename Fn>
    auto partitionner(Particles& particles, Box<T, dim> const& box, Fn&& fn)
    {
        return partitionner(particles.begin(), particles.end(), box, std::forward<Fn>(fn));
    }


    template<typename ParticleIterator, template<typename> typename Vector, typename Box>
    auto partitionner(ParticleIterator begin, ParticleIterator end, Vector<Box> const& boxes)
    {
        std::vector<BoxRange<Box, ParticleIterator>> iterators;
        auto pivot = begin;
        for (auto const& box : boxes)
            pivot = iterators.emplace_back(partitionner(pivot, end, box)).end();
        return iterators;
    }


} // namespace core
} // namespace PHARE
#endif
