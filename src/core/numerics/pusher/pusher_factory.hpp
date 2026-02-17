#ifndef PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_HPP
#define PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_HPP

#include "boris.hpp"

#include <memory>
#include <string>
#include <cstddef>


namespace PHARE
{
namespace core
{
    class PusherFactory
    {
    public:
        template<std::size_t dim, typename ParticleRange, typename Electromag,
                 typename Interpolator, typename BoundaryCondition, typename GridLayout>
        static auto makePusher(std::string pusherName)
        {
            // return type auto will fail if there's ever a second pusher

            if (pusherName == "modified_boris")
            {
                return std::make_unique<BorisPusher<dim, ParticleRange, Electromag, Interpolator,
                                                    BoundaryCondition, GridLayout>>();
            }

            throw std::runtime_error("Error : Invalid Pusher name");
        }
    };
} // namespace core
} // namespace PHARE

#endif // PUSHER_FACTORY_HPP
