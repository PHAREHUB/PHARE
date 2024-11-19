#ifndef PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_HPP
#define PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_HPP

#include <cstddef>
#include <memory>
#include <string>

#include "pusher.hpp"
#include "boris.hpp"
#include "fusher.hpp"
#include "chuck_norris.hpp"

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

namespace PHARE::core::other
{
class PusherFactory
{
public:
    template<std::size_t dim, typename ParticleRange, typename Electromag, typename Interpolator,
             typename BoundaryCondition, typename GridLayout>
    static auto makePusher(std::string pusherName)
    {
        if (pusherName == "modified_boris")
        {
            return std::make_unique<BorisPusher<dim, ParticleRange, Electromag, Interpolator,
                                                BoundaryCondition, GridLayout>>();
        }

        throw std::runtime_error("Error : Invalid Pusher name");
    }
};

} // namespace PHARE::core::other

#endif // PUSHER_FACTORY_HPP
