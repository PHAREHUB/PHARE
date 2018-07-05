#ifndef PHARE_PARTICLES_DATA_SPLIT_H
#define PHARE_PARTICLES_DATA_SPLIT_H

#include "data/particles/particles_data.h"
#include "split.h"
#include "tools/amr_utils.h"

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/pdat/CellOverlap.h>

#include <functional>

namespace PHARE
{
template<std::size_t interpOrder>
std::size_t constexpr ghostWidthForParticles()
{
    return (interpOrder % 2 == 0 ? interpOrder / 2 + 1 : (interpOrder + 1) / 2);
}

enum class ParticlesDataSplitType {
    coarseBoundary,
    interior,
    coarseBoundaryOld,
    coarseBoundaryNew
};


template<std::size_t dim, std::size_t interpOrder, ParticlesDataSplitType splitType,
         std::size_t refinedParticleNbr, typename SplitT>
class ParticlesDataSplitOperator : public SAMRAI::hier::RefineOperator
{
public:
    ParticlesDataSplitOperator()
        : SAMRAI::hier::RefineOperator{"ParticlesDataSplit_" + splitName_(splitType)}
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




    /** @brief perform a split and keep those that are inside a fineOverlap
     *
     */
    void refine(SAMRAI::hier::Patch& destination, SAMRAI::hier::Patch const& source,
                int const destinationComponent, int const sourceComponent,
                SAMRAI::hier::BoxOverlap const& fineOverlap,
                SAMRAI::hier::IntVector const& ratio) const
    {
        // For the particles we index them as a CellIndex (for the iCell)
        // ie the particles in the iCell live between lower left node of the iCell
        // and upper right of the same iCell
        auto const& destinationFieldOverlap
            = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(fineOverlap);


        // We then need to get our ParticlesData from the patch
        auto destinationParticlesData = std::dynamic_pointer_cast<ParticlesData<dim>>(
            destination.getPatchData(destinationComponent));

        auto const sourceParticlesData
            = std::dynamic_pointer_cast<ParticlesData<dim>>(source.getPatchData(sourceComponent));

        // Finnaly we need the cartesion geometry of both patch.
        auto patchGeomDestination = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            destination.getPatchGeometry());

        auto patchGeomSource = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            source.getPatchGeometry());



        TBOX_ASSERT(destinationParticlesData);
        TBOX_ASSERT(sourceParticlesData);
        TBOX_ASSERT(patchGeomDestination);
        TBOX_ASSERT(patchGeomSource);



        // We have a correct data type, we can now perform the refine
        refine_(*destinationParticlesData, *sourceParticlesData, destinationFieldOverlap, ratio,
                *patchGeomDestination, *patchGeomSource);
    }



private:
    std::string splitName_(ParticlesDataSplitType splitTypeUsed)
    {
        switch (splitTypeUsed)
        {
            case ParticlesDataSplitType::coarseBoundary: return "coarseBoundary";
            case ParticlesDataSplitType::interior: return "interior";
            case ParticlesDataSplitType::coarseBoundaryOld: return "coarseBoundaryOld";
            case ParticlesDataSplitType::coarseBoundaryNew: return "coarseBoundaryNew";
            default: throw std::runtime_error("End of enum class possible range");
        }
    }




    /** @brief given a two ParticlesData (destination and source),
     * an overlap , a ratio and the geometry of both patch, perform the
     * splitting of coarse particules onto the destination patch
     *
     */
    void refine_(ParticlesData<dim>& destinationParticlesData,
                 ParticlesData<dim> const& sourceParticlesData,
                 SAMRAI::pdat::CellOverlap const& destinationFieldOverlap,
                 SAMRAI::hier::IntVector const& ratio,
                 SAMRAI::geom::CartesianPatchGeometry const& patchGeomDest,
                 SAMRAI::geom::CartesianPatchGeometry const& patchGeomSrc) const
    {
        auto const& destinationBoxes = destinationFieldOverlap.getDestinationBoxContainer();


        // We get the reference of particles data on : interior, ghost , coarseBoundary
        // coarseBoundaryOld , coarseBoundaryNew

        auto const& sourceInteriorParticles = sourceParticlesData.domainParticles;
        auto const sourceGhostParticles     = sourceParticlesData.ghostParticles;

        auto& destinationCoarseBoundaryParticles = destinationParticlesData.coarseToFineParticles;
        auto& destinationDomainParticles         = destinationParticlesData.domainParticles;

        auto& destinationCoarseBoundaryOldParticles
            = destinationParticlesData.coarseToFineParticlesOld;
        auto& destinationCoarseBoundaryNewParticles
            = destinationParticlesData.coarseToFineParticlesNew;

        // We get the source box that contains ghost region in order to get local index later
        // same for destinationGhostBox and destinationDomainBox the later will allow to get an
        // index relative to the interior
        auto const& sourceGhostBox       = sourceParticlesData.getGhostBox();
        auto const& destinationGhostBox  = destinationParticlesData.getGhostBox();
        auto const& destinationDomainBox = destinationParticlesData.getBox();



        auto computeRatio = [&ratio]() {
            Point<int32, dim> pointRatio;

            for (auto iDir = dirX; iDir < dim; ++iDir)
            {
                pointRatio[iDir] = ratio[iDir];
            }
            return pointRatio;
        };

        SplitT split{computeRatio(), refinedParticleNbr};

        // The PatchLevelFillPattern had compute boxes that correspond to the expected filling.
        // In case of a coarseBoundary it will most likely give multiple boxes
        // in case of interior, this will be just one boxe usually
        for (auto const& destinationBox : destinationBoxes)
        {
            auto localDestinationBox = destinationBox;

            auto localSourceBox = AMRToLocal(sourceGhostBox, sourceGhostBox);

            localDestinationBox = AMRToLocal(
                static_cast<std::add_const_t<decltype(localDestinationBox)>>(localDestinationBox),
                destinationGhostBox);


            // this is same as isIn but with samrai boxes
            auto isInBox = [](auto const& particle, auto const& localSourceBox) {
                //
                auto isIn1D = [](auto const& pos, auto const& lower, auto const& upper) {
                    return pos >= lower && pos <= upper;
                };

                bool isIn{true};

                for (auto iDir = dirX; iDir < dim; ++iDir)
                {
                    isIn = isIn
                           && isIn1D(particle.iCell[iDir], localSourceBox.lower(iDir),
                                     localSourceBox.upper(iDir));
                }
                return isIn;
            };




            // Now we need to get the physical position of the boundary of the destinationBox.
            // for that we need to consider the domainBox of the destinationPatchData as the
            // reference Box

            auto* dxDest     = patchGeomDest.getDx();
            auto* xLowerDest = patchGeomDest.getXLower();

            auto* dxSrc     = patchGeomSrc.getDx();
            auto* xLowerSrc = patchGeomSrc.getXLower();


            Box<double, dim> physicalBoxDestination;

            auto destinationBoxLocalToDomain = AMRToLocal(
                static_cast<std::add_const_t<decltype(destinationBox)>>(destinationBox),
                destinationDomainBox);



            for (auto iDir = dirX; iDir < dim; ++iDir)
            {
                physicalBoxDestination.lower[iDir]
                    = xLowerDest[iDir] + dxDest[iDir] * destinationBoxLocalToDomain.lower(iDir);

                if constexpr (splitType == ParticlesDataSplitType::interior)
                {
                    if constexpr (interpOrder == 1)
                    {
                        physicalBoxDestination.upper[iDir]
                            = xLowerDest[iDir]
                              + dxDest[iDir]
                                    * (destinationBoxLocalToDomain.upper(iDir) + 1
                                       + ghostWidthForParticles<interpOrder>());
                    }
                    else if constexpr (interpOrder == 2 || interpOrder == 3)
                    {
                        physicalBoxDestination.upper[iDir]
                            = xLowerDest[iDir]
                              + dxDest[iDir]
                                    * (destinationBoxLocalToDomain.upper(iDir) + 1
                                       + ghostWidthForParticles<interpOrder>() + 1);
                    }
                }
                else
                {
                    physicalBoxDestination.upper[iDir]
                        = xLowerDest[iDir]
                          + dxDest[iDir] * (destinationBoxLocalToDomain.upper(iDir) + 1);
                }
            }


            Point<double, dim> originDest;
            Point<double, dim> originSrc;

            for (auto iDir = dirX; iDir < dim; ++iDir)
            {
                originDest[iDir] = xLowerDest[iDir]
                                   - sourceParticlesData.getGhostCellWidth()[iDir] * dxDest[iDir];
                originSrc[iDir]
                    = xLowerSrc[iDir] - sourceParticlesData.getGhostCellWidth()[iDir] * dxSrc[iDir];
            }


            // this will set coarse particle position as a
            // refined particles position (by adapting icell and delta)
            auto particleAtRefinedPosition
                = [&sourceGhostBox, &destinationGhostBox, &ratio, &originDest, &originSrc, dxDest,
                   dxSrc](auto& particle) {
                      //

                      Point<double, dim> position{positionAsPoint(particle, dxSrc, originSrc)};

                      for (auto iDir = dirX; iDir < dim; ++iDir)
                      {
                          double normalizedPos = (position[iDir] - originDest[iDir]) / dxDest[iDir];

                          particle.iCell[iDir] = static_cast<int>(normalizedPos);
                          particle.delta[iDir]
                              = normalizedPos - static_cast<double>(particle.iCell[iDir]);
                      }
                  };




            auto isInSplit = [&physicalBoxDestination, dxDest, &originDest](auto const& particle) {
                Point<double, dim> particlesPosition{positionAsPoint(particle, dxDest, originDest)};
                return isIn(particlesPosition, physicalBoxDestination);
            };

            // Since we are in a temporary space, we may have to copy information
            // from ghost region as well. This operator will perform the split
            // on particles in domain and ghost zone, and put the split particles
            // in the coarseToFineParticles.
            std::array<std::remove_reference_t<decltype(sourceInteriorParticles)>*, 2>
                particlesArrays{{&sourceInteriorParticles, &sourceGhostParticles}};



            // We loop over interiorParticles and ghostParticles
            // for each particles, if they are in the localSourceBox
            // then we shift them in the destination space (coarseAtRefinedPosition)
            // case 1-3
            // if we want to split on the coarseBoundary we check if the coarseParticle is
            // a candidate to split. If it is the case we split in a temporary vector
            // then for each particles in the temporary vector, we check if they are inside
            // the desired region, and put them in.
            // case 4
            // if we want to split on the interior, we split each particles in a temporary vector
            // and keep the ones that enter the interior domains.


            for (auto const& sourceParticlesArray : particlesArrays)
            {
                for (auto const& particle : *sourceParticlesArray)
                {
                    //
                    auto particlePosition{positionAsPoint(particle, dxSrc, originSrc)};

                    if (isInBox(particle, localSourceBox))
                    {
                        //
                        auto particleRefinedPos = particle;
                        particleAtRefinedPosition(particleRefinedPos);

                        bool constexpr isCoarseBoundarySplitType
                            = splitType == ParticlesDataSplitType::coarseBoundary
                              || splitType == ParticlesDataSplitType::coarseBoundaryOld
                              || splitType == ParticlesDataSplitType::coarseBoundaryNew;




                        if constexpr (isCoarseBoundarySplitType)
                        {
                            if (this->isCandidateForSplit_(particleRefinedPos,
                                                           physicalBoxDestination, dxDest,
                                                           originDest, dxSrc))
                            {
                                std::vector<Particle<dim>> splittedParticles;

                                auto particlePosition{
                                    positionAsPoint(particleRefinedPos, dxDest, originDest)};

                                split(particleRefinedPos, splittedParticles);

                                if constexpr (splitType == ParticlesDataSplitType::coarseBoundary)
                                {
                                    std::copy_if(
                                        std::begin(splittedParticles), std::end(splittedParticles),
                                        std::back_inserter(destinationCoarseBoundaryParticles),
                                        isInSplit);
                                }
                                else if constexpr (splitType
                                                   == ParticlesDataSplitType::coarseBoundaryOld)
                                {
                                    //
                                    std::copy_if(
                                        std::begin(splittedParticles), std::end(splittedParticles),
                                        std::back_inserter(destinationCoarseBoundaryOldParticles),
                                        isInSplit);
                                }
                                else //  splitType is coarseBoundaryNew
                                {
                                    //
                                    std::copy_if(
                                        std::begin(splittedParticles), std::end(splittedParticles),
                                        std::back_inserter(destinationCoarseBoundaryNewParticles),
                                        isInSplit);
                                }
                            }
                        }
                        else // splitType == ParticlesDataSplitType::interior
                        {
                            std::vector<Particle<dim>> splittedParticles;
                            split(particleRefinedPos, splittedParticles);

                            std::copy_if(std::begin(splittedParticles), std::end(splittedParticles),
                                         std::back_inserter(destinationDomainParticles), isInSplit);
                        }
                    } // is in source selected box
                }     // loop on particle
            }         // loop on particelesArrays
        }             // loop on destination box
    }




    template<typename Particle>
    bool isCandidateForSplit_(Particle const& particle,
                              Box<double, dim> const& physicalBoxDestination, double const* dxDest,
                              Point<double, dim> const& originDest, double const* dxSrc) const
    {
        // Given the physicalBoxDestination where splitted particles must be keep
        // return true if a coarseParticles may split in the region
        Point<double, dim> maxDistance;

        Point<double, dim> particlePosition{positionAsPoint(particle, dxDest, originDest)};
        Point<double, dim> distanceFromLower;
        Point<double, dim> distanceFromUpper;

        Box<double, dim> candidateBoxForSplit;

        // the maximum distance is ((interpOrder ) / 2. ) * dxCoarse

        for (auto iDir = dirX; iDir < dim; ++iDir)
        {
            maxDistance[iDir] = (interpOrder / 2.) * dxSrc[iDir];
        }


        for (auto iDir = dirX; iDir < dim; ++iDir)
        {
            distanceFromLower[iDir]
                = std::abs(particlePosition[iDir] - physicalBoxDestination.lower[iDir]);
            distanceFromUpper[iDir]
                = std::abs(particlePosition[iDir] - physicalBoxDestination.upper[iDir]);
        }

        bool mustSplit = true;

        for (auto iDir = dirX; iDir < dim; ++iDir)
        {
            mustSplit = mustSplit
                        && (distanceFromLower[iDir] <= maxDistance[iDir]
                            || distanceFromUpper[iDir] <= maxDistance[iDir]);
        }
        return mustSplit;
    }
};


} // namespace PHARE

#endif
