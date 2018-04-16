#ifndef PHARE_CORE_UTILITIES_PARTICLE_SELECTOR_PARTICLE_SELECTOR_H
#define PHARE_CORE_UTILITIES_PARTICLE_SELECTOR_PARTICLE_SELECTOR_H

#include <cstddef>

#include "data/particles/particle.h"
#include "utilities/box/box.h"

namespace PHARE
{
template<typename BoxCollection>
class ParticleSelector
{
public:
    explicit ParticleSelector(BoxCollection boxes)
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
    return ParticleSelector<BoxContainer>{std::forward<BoxContainer>(boxes)};
}




} // namespace PHARE


#endif
