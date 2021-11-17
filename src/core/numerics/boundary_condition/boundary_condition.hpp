#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_BOUNDARY_CONDITION_HPP

#include <cstddef>

#include "core/utilities/box/box.hpp"
#include "core/utilities/partitionner/partitionner.hpp"


namespace PHARE
{
namespace core
{
    template<std::size_t dim, std::size_t interpOrder>
    class BoundaryCondition
    {
    public:
        void setBoundaryBoxes(std::vector<Box<int, dim>> boxes)
        {
            boundaryBoxes_ = std::move(boxes);
        }

        template<typename ParticleIterator>
        ParticleIterator applyOutgoingParticleBC(ParticleIterator begin, ParticleIterator end)
        {
            // TODO
            return ParticleIterator{begin}; // place holder
        }


    private:
        std::vector<Box<int, dim>> boundaryBoxes_;
    };

} // namespace core
} // namespace PHARE



#endif
