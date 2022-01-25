#ifndef PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_HPP
#define PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_HPP

#include <cstddef>
#include <memory>
#include <string>

#include "boris.hpp"
#include "pusher.hpp"

namespace PHARE
{
namespace core
{
    class PusherFactory
    {
    public:
        template<std::size_t dim, typename ParticleIterator, typename Particle_t,
                 typename Electromag, typename Interpolator, typename BoundaryCondition,
                 typename GridLayout>
        static auto makePusher(std::string pusherName)
        {
            if (pusherName == "modified_boris")
            {
                return std::make_unique<BorisPusher<dim, ParticleIterator, Particle_t, Electromag,
                                                    Interpolator, BoundaryCondition, GridLayout>>();
            }

            throw std::runtime_error("Error : Invalid Pusher name");
        }
    };
} // namespace core
} // namespace PHARE

#endif // PUSHER_FACTORY_HPP
