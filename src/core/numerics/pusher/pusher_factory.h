#ifndef PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_H
#define PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_H

#include <cstddef>
#include <string>
#include <memory>

#include "pusher.h"
#include "boris.h"

namespace PHARE
{

class PusherFactory
{
public:
    template<std::size_t dim, typename ParticleIterator, typename Electromag, typename Interpolator,
             typename ParticleSelector, typename BoundaryCondition>
    static auto makePusher(std::string pusherName)
    {
        if (pusherName == "boris")
        {
            return std::make_unique<BorisPusher<dim, ParticleIterator, Electromag, Interpolator,
                    ParticleSelector, BoundaryCondition>>();
        }

        else
        {
            std::runtime_error("Error : Invalid Pusher name");
        }
    }
};

}

#endif // PUSHER_FACTORY_H
