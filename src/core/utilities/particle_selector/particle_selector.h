#ifndef PHARE_CORE_UTILITIES_PARTICLE_SELECTOR_PARTICLE_SELECTOR_H
#define PHARE_CORE_UTILITIES_PARTICLE_SELECTOR_PARTICLE_SELECTOR_H

#include <cstddef>

#include "core/data/particles/particle.h"
#include "core/utilities/box/box.h"

namespace PHARE
{
namespace core
{
    template<typename BoxCollection>
    class ParticleSelectorT
    {
    public:
        explicit ParticleSelectorT(BoxCollection boxes)
            : boxes_{boxes}
        {
        }

        template<typename Particle>
        bool operator()(Particle const& particle) const
        {
            auto point = cellAsPoint(particle);
            return isIn(point, boxes_);
        }


    private:
        BoxCollection boxes_;
    };


    template<typename BoxContainer>
    auto makeSelector(BoxContainer&& boxes)
    {
        return ParticleSelectorT<BoxContainer>{std::forward<BoxContainer>(boxes)};
    }

} // namespace core


} // namespace PHARE


#endif
