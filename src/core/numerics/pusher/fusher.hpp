#ifndef PHARE_CORE_NUMERICS_PUSHER_FUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_FUSHER_HPP

#include <cstddef>
#include <utility>
#include <functional>
#include <type_traits>

// #include "core/utilities/range/range.hpp"
// #include "core/data/particles/particle.hpp"

namespace PHARE::core::other
{
template<std::size_t dim, typename ParticleArray_t, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout, bool UseParticleIndex = true>
class Pusher
{
protected:
    static auto constexpr dimension = GridLayout::dimension;
    using Particle_t                = ParticleArray_t::value_type;
    using OnChangeIcell             = std::function<bool(ParticleArray_t&, Particle_t&,
                                             std::array<int, dim> const&, std::size_t const)>;

public:
    virtual void move(ParticleArray_t& array, Electromag const& emFields, double mass,
                      Interpolator& interpolator, GridLayout const& layout,
                      OnChangeIcell&& icellchanger, bool ghosts = false)
        = 0;


    virtual void setMeshAndTimeStep(std::array<double, dim> ms, double ts) = 0;

    virtual ~Pusher() {}
};

} // namespace PHARE::core::other

#endif
