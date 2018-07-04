#ifndef PHARE_PARTICLES_DATA_SPLIT_ON_COARSE_BOUNDARY_H
#define PHARE_PARTICLES_DATA_SPLIT_ON_COARSE_BOUNDARY_H

#include "data/particles/particles_data.h"
#include "tools/amr_utils.h"

#include <SAMRAI/geom/CartesianPatchGeometry.h>
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

enum class ParticlesDataSplitType { coarseBoundary, interior, coarseBoundary1, coarseBoundary2 };

std::string inline splitName(ParticlesDataSplitType splitType)
{
    switch (splitType)
    {
        case ParticlesDataSplitType::coarseBoundary: return "coarseBoundary";
        case ParticlesDataSplitType::interior: return "interior";
        case ParticlesDataSplitType::coarseBoundary1: return "coarseBoundary1";
        case ParticlesDataSplitType::coarseBoundary2: return "coarseBoundary2";
        default: throw std::runtime_error("End of enum class possible range");
    }
}

template<std::size_t dim, std::size_t interpOrder, ParticlesDataSplitType splitType>
class ParticlesDataSplitOperator : public SAMRAI::hier::RefineOperator
{
public:
    ParticlesDataSplitOperator()
        : SAMRAI::hier::RefineOperator{"ParticlesDataSplit_" + splitName(splitType)}
    {
    }

    virtual ~ParticlesDataSplitOperator() = default;

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

        auto const& sourceGhostBox       = sourceParticlesData.getGhostBox();
        auto const& destinationGhostBox  = destinationParticlesData.getGhostBox();
        auto const& destinationDomainBox = destinationParticlesData.getBox();

        for (auto const& destinationBox : destinationBoxes)
        {
            auto sourceBox = destinationBox;

            sourceBox.coarsen(ratio);

            auto localDestinationBox = destinationBox;

            if constexpr (splitType == ParticlesDataSplitType::coarseBoundary)
            {
                auto growthVector = SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dim},
                                                            ghostWidthForParticles<interpOrder>()};

                sourceBox.grow(growthVector);
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
                particle.iCell = localToAMR(particle.iCell, sourceGhostBox);
                particle.iCell = refinedPosition(particle.iCell, ratio);
                particle.iCell = AMRToLocal(particle.iCell, destinationGhostBox);


                auto shiftDeltaAndiCell = [&particle](auto direction) {
                    // TODO find better
                    constexpr double half = 0.5;
                    if (particle.delta[direction] >= half)
                    {
                        ++particle.iCell[direction];
                        particle.delta[direction] -= half;
                    }
                };

                shiftDeltaAndiCell(dirX);
                if constexpr (dim > 1)
                {
                    shiftDeltaAndiCell(dirY);
                }
                if constexpr (dim > 2)
                {
                    shiftDeltaAndiCell(dirZ);
                }
            };

            auto pGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
                destination.getPatchGeometry());

            auto* dx     = pGeom->getDx();
            auto* xLower = pGeom->getXLower();

            Point<double, dim> physicalLowerDestination;
            Point<double, dim> physicalUpperDestination;

            auto destinationBoxLocalToDomain = AMRToLocal(
                static_cast<std::add_const_t<decltype(destinationBox)>>(destinationBox),
                destinationDomainBox);

            physicalLowerDestination[dirX]
                = xLower[dirX] + dx[dirX] * destinationBoxLocalToDomain.lower(dirX);
            physicalUpperDestination[dirX]
                = xLower[dirX] + dx[dirX] * (destinationBoxLocalToDomain.upper(dirX) + 1);

            if constexpr (dim > 1)
            {
                physicalLowerDestination[dirY]
                    = xLower[dirY] + dx[dirY] * destinationBoxLocalToDomain.lower(dirY);

                physicalUpperDestination[dirY]
                    = xLower[dirY] + dx[dirY] * (destinationBoxLocalToDomain.upper(dirY) + 1);
            }
            if constexpr (dim > 2)
            {
                physicalLowerDestination[dirZ]
                    = xLower[dirZ] + dx[dirZ] * destinationBoxLocalToDomain.lower(dirZ);

                physicalUpperDestination[dirZ]
                    = xLower[dirZ] + dx[dirZ] * (destinationBoxLocalToDomain.upper(dirZ) + 1);
            }

            [[maybe_unused]] auto isCandidateForSplit = [&physicalLowerDestination,
                                                         &physicalUpperDestination, dx,
                                                         xLower](auto const& particle) {
                if constexpr (dim == 1)
                {
                    if constexpr (interpOrder == 1)
                    {
                        double maxDistanceX       = 0.5 * dx[dirX];
                        double particlesPositionX = xLower[dirX] + particle.iCell[dirX] * dx[dirX]
                                                    + particle.delta[dirX] * dx[dirX];
                        double distanceFromLowerX
                            = std::abs(particlesPositionX - physicalLowerDestination[dirX]);
                        double distanceFromUpperX
                            = std::abs(particlesPositionX - physicalUpperDestination[dirX]);

                        return (distanceFromLowerX <= maxDistanceX
                                || distanceFromUpperX <= maxDistanceX);
                    }
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

                        if constexpr (splitType == ParticlesDataSplitType::coarseBoundary)
                        {
                            if (isCandidateForSplit(shiftedParticle))
                            {
                                destinationCoarseBoundaryParticles.push_back(shiftedParticle);
                            }
                        }
                        else if constexpr (splitType == ParticlesDataSplitType::interior)
                        {
                            if (isInBox(shiftedParticle, localDestinationBox))
                            {
                                destinationDomainParticles.push_back(shiftedParticle);
                            }
                        }
                    } // is in source selected box
                }     // loop on particle
            }         // loop on particelesArrays
        }             // loop on destination box
    }
};


} // namespace PHARE

#endif
