#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_BOUNDARY_CONDITION_H
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_BOUNDARY_CONDITION_H

#include <cstddef>

#include "core/utilities/box/box.h"
#include "core/utilities/partitionner/partitionner.h"


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
            // TODO loop while last partition not equel to begin.
            auto partitions = partitionner(begin, end, boundaryBoxes_);

            // applyBC for each box.
            // end while loop

            // particles to delete are between te last iterator
            // of 'partitions' and the end.
            return partitions[boundaryBoxes_.size()];
        }


    private:
        std::vector<Box<int, dim>> boundaryBoxes_;
    };

} // namespace core
} // namespace PHARE



#endif
