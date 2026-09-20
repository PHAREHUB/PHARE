#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PARTITIONER_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PARTITIONER_HPP

#include "core/def.hpp"
#include "core/data/vector.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/partitionner/partitionner.hpp"
#include "core/utilities/range/range.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace PHARE::core
{

// --- public API ---

namespace detail
{
    // hides the per layout/alloc impl behind operator()/notIn() - see implementation section below
    template<auto layout_mde, auto alloc_mde, typename ParticleArray>
    struct ParticleArrayPartitioner;
} // namespace detail

template<typename ParticleArray>
struct ParticleArrayPartitioner
{
    using box_t = Box<int, ParticleArray::dimension>;
    using _impl = detail::ParticleArrayPartitioner<ParticleArray::layout_mode,
                                                   ParticleArray::alloc_mode, ParticleArray>;

    ParticleArrayPartitioner(ParticleArray& particles_)
        : ParticleArrayPartitioner{particles_, 0, particles_.size()}
    {
    }
    ParticleArrayPartitioner(ParticleArray& particles_, std::size_t start_, std::size_t end_)
        : particles{particles_}
        , start{start_}
        , end{end_}
    {
    }

    template<typename... Args>
    auto operator()(Args&&... args)
    {
        return _impl{particles, start, end}(std::forward<Args>(args)...);
    }

    template<typename... Args>
    auto notIn(Args&&... args)
    {
        return _impl{particles, start, end}.notIn(std::forward<Args>(args)...);
    }

    ParticleArray& particles;
    std::size_t start, end;
};


template<typename Particles_t>
auto partition_particles(Particles_t& particles, auto const box)
{
    return ParticleArrayPartitioner<Particles_t>{particles}(box);
}
template<typename Particles_t>
auto partition_particles_not_in(Particles_t& particles, auto const box)
{
    return ParticleArrayPartitioner<Particles_t>{particles}.notIn(box);
}

} // namespace PHARE::core


// --- implementations ---

namespace PHARE::core::detail
{

template<auto layout_mde, auto alloc_mde, typename ParticleArray>
struct ParticleArrayPartitioner
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));
    static_assert(dependent_false_v<ParticleArray>,
                  "ParticleArrayPartitioner not implemented for this layout/alloc permutation");
};


using enum LayoutMode;
using enum AllocatorMode;

// shared by the layouts that partition via the generic iterator-range partitionner():
// AoS, AoSPCTS, AoSPC
template<typename ParticleArray>
struct GenericParticleArrayPartitioner
{
    using box_t    = Box<int, ParticleArray::dimension>;
    using iterator = decltype(std::declval<ParticleArray>().begin());
    using range_t  = BoxRange<box_t, iterator>;

    GenericParticleArrayPartitioner(ParticleArray& array_)
        : GenericParticleArrayPartitioner{array_, 0, array_.size()}
    {
    }
    GenericParticleArrayPartitioner(ParticleArray& array_, std::size_t start_, std::size_t end_)
        : array{array_}
        , start{start_}
        , end{end_}
    {
    }

    auto operator()(box_t const& box)
    {
        return partitionner(array.begin() + start, array.begin() + end, box);
    }

    auto notIn(box_t const& box)
    {
        return partitionner(array.begin() + start, array.begin() + end, box,
                            [=](auto& part) { return !isIn(cellAsPoint(part.iCell()), box); });
    }

    template<std::size_t S>
    auto operator()(std::array<box_t, S> const& boxes)
    {
        static_assert(S > 0);
        auto iterators = generate_from(
            [&](auto const& box) { return range_t{box, array.begin(), array.begin()}; }, boxes);
        for (std::size_t i = 0; i < boxes.size(); ++i)
            start += (iterators[i] = (*this)(boxes[i])).size();
        return iterators;
    }

    auto operator()(std::vector<box_t> const& boxes)
    {
        std::vector<range_t> iterators;
        for (auto const& box : boxes)
            start += iterators.emplace_back((*this)(box)).size();
        return iterators;
    }

    ParticleArray& array;
    std::size_t start, end;
};


template<typename ParticleArray>
struct ParticleArrayPartitioner<AoS, CPU, ParticleArray>
    : GenericParticleArrayPartitioner<ParticleArray>
{
    using GenericParticleArrayPartitioner<ParticleArray>::GenericParticleArrayPartitioner;
};

template<typename ParticleArray>
struct ParticleArrayPartitioner<AoSPCTS, CPU, ParticleArray>
    : GenericParticleArrayPartitioner<ParticleArray>
{
    using GenericParticleArrayPartitioner<ParticleArray>::GenericParticleArrayPartitioner;
};

template<typename ParticleArray>
struct ParticleArrayPartitioner<AoSPC, CPU, ParticleArray>
    : GenericParticleArrayPartitioner<ParticleArray>
{
    using GenericParticleArrayPartitioner<ParticleArray>::GenericParticleArrayPartitioner;
};


template<typename ParticleArray>
struct ParticleArrayPartitioner<AoSMapped, CPU, ParticleArray>
{
    using box_t    = Box<int, ParticleArray::dimension>;
    using iterator = decltype(std::declval<ParticleArray>().begin());
    using range_t  = BoxRange<box_t, iterator>;

    auto operator()(box_t const& box)
    {
        return array.partition([&](auto& part) { return isIn(cellAsPoint(part.iCell()), box); });
    }

    auto notIn(box_t const& box)
    {
        return array.partition([&](auto& part) { return !isIn(cellAsPoint(part.iCell()), box); });
    }

    template<std::size_t S>
    auto operator()(std::array<box_t, S> const& boxes)
    {
        static_assert(S > 0);
        auto iterators = generate_from(
            [&](auto const& box) { return range_t{box, array.begin(), array.begin()}; }, boxes);
        for (std::size_t i = 0; i < boxes.size(); ++i)
            start += (iterators[i] = (*this)(boxes[i])).size();
        return iterators;
    }

    auto operator()(std::vector<box_t> const& boxes)
    {
        std::vector<range_t> iterators;
        for (auto const& box : boxes)
            start += iterators.emplace_back((*this)(box)).size();
        return iterators;
    }

    ParticleArray& array;
    std::size_t start = 0, end = array.size();
};


// SoA is internal-only plumbing (HDF5/restart write, python interop) and is never partitioned -
// kept for symmetry, same reasoning as particle_array_sorter.hpp
template<typename ParticleArray>
struct ParticleArrayPartitioner<SoA, CPU, ParticleArray>
{
    using box_t    = Box<int, ParticleArray::dimension>;
    using iterator = decltype(std::declval<ParticleArray>().begin());
    using range_t  = BoxRange<box_t, iterator>;

    range_t operator()(box_t const&) { throw std::runtime_error("no partition for SoA"); }
    range_t notIn(box_t const&) { throw std::runtime_error("no partition for SoA"); }

    template<std::size_t S>
    auto operator()(std::array<box_t, S> const& boxes)
    {
        static_assert(S > 0);
        auto iterators = generate_from(
            [&](auto const& box) { return range_t{box, array.begin(), array.begin()}; }, boxes);
        for (std::size_t i = 0; i < boxes.size(); ++i)
            start += (iterators[i] = (*this)(boxes[i])).size();
        return iterators;
    }

    auto operator()(std::vector<box_t> const& boxes)
    {
        std::vector<range_t> iterators;
        for (auto const& box : boxes)
            start += iterators.emplace_back((*this)(box)).size();
        return iterators;
    }

    ParticleArray& array;
    std::size_t start = 0, end = array.size();
};


} // namespace PHARE::core::detail

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PARTITIONER_HPP */
