#ifndef PHARE_PARTICLES_DATA_SPLIT_ON_COARSE_BOUNDARY_H
#define PHARE_PARTICLES_DATA_SPLIT_ON_COARSE_BOUNDARY_H

#include "data/particles/particles_data.h"
#include "tools/amr_utils.h"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/pdat/CellOverlap.h>

namespace PHARE
{
template<std::size_t interpOrder>
std::size_t constexpr ghostWidthForParticles()
{
    return interpOrder % 2 == 0 ? interpOrder / 2 + 1 : (interpOrder + 1) / 2;
}

template<std::size_t dim, std::size_t interpOrder>
class ParticlesDataSplitOnCoarseBoundary : public SAMRAI::hier::RefineOperator
{
public:
    explicit ParticlesDataSplitOnCoarseBoundary(bool refineOnBorderOnly)
        : SAMRAI::hier::RefineOperator{"ParticlesDataSplitOnCoarseBoundary_"
                                       + std::to_string(refineOnBorderOnly)}
        , refineOnBorderOnly_{refineOnBorderOnly}
    {
    }

    virtual ~ParticlesDataSplitOnCoarseBoundary() = default;

    /** @brief a priority of 0 means that this operator
     * will be applied first
     */
    int getOperatorPriority() const { return 0; }

    SAMRAI::hier::IntVector getStencilWidth(SAMRAI::tbox::Dimension const& dimension) const
    {
        return SAMRAI::hier::IntVector{dimension, ghostWidthForParticles<interpOrder>()};
    }

    void refine(SAMRAI::hier::Patch& destination, SAMRAI::hier::Patch const& source,
                int const destinationComponent, int const sourceComponent,
                SAMRAI::hier::BoxOverlap const& fineOverlap,
                SAMRAI::hier::IntVector const& ratio) const
    {
        auto const& destinationFieldOverlap
            = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(fineOverlap);

        auto const& destinationBoxes = destinationFieldOverlap.getDestinationBoxContainer();

        auto& destinationParticlesData = *std::dynamic_pointer_cast<ParticlesData<dim>>(
            destination.getPatchData(destinationComponent));

        auto const& sourceParticlesData
            = *std::dynamic_pointer_cast<ParticlesData<dim>>(source.getPatchData(sourceComponent));


        auto const& sourceGeometry = source.getPatchGeometry();


        TBOX_ASSERT(destinationFieldData);
        TBOX_ASSERT(sourceFieldData);

        auto const& sourceInteriorParticles = sourceParticlesData.domainParticles;
        auto const sourceGhostParticles     = sourceParticlesData.ghostParticles;

        auto& destinationCoarseBoundaryParticles = destinationParticlesData.coarseToFineParticles;
        auto& destinationDomainParticles         = destinationParticlesData.domainParticles;

        auto const& sourceGhostBox      = sourceParticlesData.getGhostBox();
        auto const& destinationGhostBox = destinationParticlesData.getGhostBox();

        for (auto const& destinationBox : destinationBoxes)
        {
            auto sourceBox = destinationBox;

            sourceBox.coarsen(ratio);

            auto localDestinationBox = destinationBox;

            // TODO if constexpr in the future
            if (refineOnBorderOnly_)
            {
                auto growthVector = SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dim},
                                                            ghostWidthForParticles<interpOrder>()};

                sourceBox.grow(growthVector);
                localDestinationBox.grow(growthVector);
            }

            sourceBox = sourceBox * sourceGhostBox;

            auto localSourceBox = AMRToLocal(
                static_cast<std::add_const_t<decltype(sourceBox)>>(sourceBox), sourceGhostBox);

            localDestinationBox = AMRToLocal(
                static_cast<std::add_const_t<decltype(localDestinationBox)>>(localDestinationBox),
                destinationGhostBox);

            auto isInBox = [](auto const& particle, auto const& localSourceBox) {
                if constexpr (dim == 1)
                {
                    return particle.iCell[dirX] >= localSourceBox.lower(dirX)
                           && particle.iCell[dirX] <= localSourceBox.upper(dirX);
                }
                else if constexpr (dim == 2)
                {
                    auto inDirX = particle.iCell[dirX] >= localSourceBox.lower(dirX)
                                  && particle.iCell[dirX] <= localSourceBox.upper(dirX);
                    auto inDirY = particle.iCell[dirY] >= localSourceBox.lower(dirY)
                                  && particle.iCell[dirY] <= localSourceBox.upper(dirY);

                    return inDirX && inDirY;
                }
                else if constexpr (dim == 3)
                {
                    auto inDirX = particle.iCell[dirX] >= localSourceBox.lower(dirX)
                                  && particle.iCell[dirX] <= localSourceBox.upper(dirX);

                    auto inDirY = particle.iCell[dirY] >= localSourceBox.lower(dirY)
                                  && particle.iCell[dirY] <= localSourceBox.upper(dirY);

                    auto inDirZ = particle.iCell[dirZ] >= localSourceBox.lower(dirZ)
                                  && particle.iCell[dirZ] <= localSourceBox.upper(dirZ);

                    return inDirX && inDirY && inDirZ;
                }
            };

            auto shiftParticle = [&sourceGhostBox, &destinationGhostBox, &ratio](auto& particle) {
                auto const& lowerSourceGhostBox      = sourceGhostBox.lower();
                auto const& lowerDestinationGhostBox = destinationGhostBox.lower();

                auto shiftDirection = [&lowerSourceGhostBox, &ratio, &lowerDestinationGhostBox,
                                       &particle](auto direction) {
                    //
                    particle.iCell[direction] += lowerSourceGhostBox(direction);
                    particle.iCell[direction] *= ratio(direction);
                    particle.iCell[direction] -= lowerDestinationGhostBox(direction);

                    // TODO: this  will be replace with a split algorithm
                    particle.iCell[direction]
                        += static_cast<int>(particle.delta[direction] + 0.5) / 2;
                };
                shiftDirection(dirX);
                if constexpr (dim > 1)
                {
                    shiftDirection(dirY);
                }
                if constexpr (dim > 2)
                {
                    shiftDirection(dirZ);
                }
            };

            // Since we are in a temporary space, we may have to copy information
            // from ghost region as well. This operator will perform the split
            // on particles in domain and ghost zone, and put the split particles
            // in the coarseToFineParticles.
            std::array<std::remove_reference_t<decltype(sourceInteriorParticles)>*, 2>
                particlesArrays{{&sourceInteriorParticles, &sourceGhostParticles}};

            for (auto const& sourceParticlesArray : particlesArrays)
            {
                for (auto const& particle : *sourceParticlesArray)
                {
                    //

                    if (isInBox(particle, localSourceBox))
                    {
                        //
                        auto shiftedParticle = particle;
                        shiftParticle(shiftedParticle);
                        if (isInBox(shiftedParticle, localDestinationBox))
                        {
                            // TODO replace with constexpr
                            if (refineOnBorderOnly_)
                            {
                                destinationCoarseBoundaryParticles.push_back(shiftedParticle);
                            }
                            else
                            {
                                destinationDomainParticles.push_back(shiftedParticle);
                            }
                        }
                    }
                }
            }

            //
        }
    }

private:
    const bool refineOnBorderOnly_;
};


} // namespace PHARE

#endif
