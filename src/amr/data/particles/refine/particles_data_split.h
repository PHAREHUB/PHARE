#ifndef PHARE_PARTICLES_DATA_SPLIT_H
#define PHARE_PARTICLES_DATA_SPLIT_H

#include "amr/data/particles/particles_data.h"
#include "amr/resources_manager/amr_utils.h"
#include "split.h"
#include "core/utilities/constants.h"

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/pdat/CellOverlap.h>

#include <functional>

namespace PHARE
{
namespace amr
{
    using core::int32;
    using core::uint32;

    enum class ParticlesDataSplitType {
        coarseBoundary,
        interior,
        coarseBoundaryOld,
        coarseBoundaryNew
    };

    template<std::size_t dim, std::size_t interpOrder, ParticlesDataSplitType splitType,
             std::size_t refinedParticleNbr, typename SplitT>
    class ParticlesRefineOperator : public SAMRAI::hier::RefineOperator
    {
    public:
        ParticlesRefineOperator()
            : SAMRAI::hier::RefineOperator{"ParticlesDataSplit_" + splitName_(splitType)}
        {
        }

        virtual ~ParticlesRefineOperator() = default;

        /** @brief a priority of 0 means that this operator
         * will be applied first
         */
        virtual int getOperatorPriority() const override { return 0; }

        virtual SAMRAI::hier::IntVector
        getStencilWidth(SAMRAI::tbox::Dimension const& dimension) const override
        {
            return SAMRAI::hier::IntVector{dimension, ghostWidthForParticles<interpOrder>()};
        }

        /** @brief perform a split and keep those that are inside a fineOverlap
         *
         */
        virtual void refine(SAMRAI::hier::Patch& destination, SAMRAI::hier::Patch const& source,
                            int const destinationComponent, int const sourceComponent,
                            SAMRAI::hier::BoxOverlap const& fineOverlap,
                            SAMRAI::hier::IntVector const& ratio) const override
        {
            // For the particles we index them as a CellIndex (for the iCell)
            // ie the particles in the iCell live between lower left node of the iCell
            // and upper right of the same iCell
            auto const& destinationFieldOverlap
                = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(fineOverlap);


            // We then need to get our ParticlesData from the patch
            auto destinationParticlesData = std::dynamic_pointer_cast<ParticlesData<dim>>(
                destination.getPatchData(destinationComponent));

            auto const sourceParticlesData = std::dynamic_pointer_cast<ParticlesData<dim>>(
                source.getPatchData(sourceComponent));

            // Finnaly we need the cartesion geometry of both patch.
            auto patchGeomDestination
                = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
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


        /** @brief given two ParticlesData (destination and source),
         * an overlap , a ratio and the geometry of both patches, perform the
         * splitting of coarse particules onto the destination patch
         */
        void refine_(ParticlesData<dim>& destParticlesData,
                     ParticlesData<dim> const& srcParticlesData,
                     SAMRAI::pdat::CellOverlap const& destFieldOverlap,
                     SAMRAI::hier::IntVector const& ratio,
                     SAMRAI::geom::CartesianPatchGeometry const& /*patchGeomDest*/,
                     SAMRAI::geom::CartesianPatchGeometry const& /*patchGeomSrc*/) const
        {
            // the source PatchData is a possible restriction of a "real" patchdata
            // so that it is the closest from the destination boxes
            // if all particles from the original source patchdata are in "domainParticles"
            // they can now be found in either domain of ghost particle arrays of this
            // temporary restriction "source" patchData
            // therefore we need references to the domain and ghost particle arrays
            auto const& srcInteriorParticles = srcParticlesData.domainParticles;
            auto const& srcGhostParticles    = srcParticlesData.patchGhostParticles;

            // the particle refine operator's job is to fill either domain (during initialization of
            // new patches) or coarse to fine boundaries (during advance), so we need references to
            // these arrays on the destination. We don't fill ghosts with this operator, they are
            // filled from exchanging with neighbor patches.
            auto const& destBoxes                = destFieldOverlap.getDestinationBoxContainer();
            auto& destCoarseBoundaryParticles    = destParticlesData.levelGhostParticles;
            auto& destDomainParticles            = destParticlesData.domainParticles;
            auto& destCoarseBoundaryOldParticles = destParticlesData.levelGhostParticlesOld;
            auto& destCoarseBoundaryNewParticles = destParticlesData.levelGhostParticlesNew;


            // We get the source box that contains ghost region in order to get local index later
            // same for destinationGhostBox and destinationDomainBox the later will allow to get an
            // index relative to the interior


            // TODO refineParticleNbr should not be runtime and SplitT should be created only once.
            // SplitT split{computeRatio(), refinedParticleNbr};

            SplitT split{core::Point<int32, dim>{ratio}, refinedParticleNbr};


            // The PatchLevelFillPattern had compute boxes that correspond to the expected filling.
            // In case of a coarseBoundary it will most likely give multiple boxes
            // in case of interior, this will be just one boxe usually
            for (auto const& destinationBox : destBoxes)
            {
                std::array<std::remove_reference_t<decltype(srcInteriorParticles)>*, 2>
                    particlesArrays{{&srcInteriorParticles, &srcGhostParticles}};


                auto isInDest = [&destinationBox](auto const& particle) //
                { return isInBox(destinationBox, particle); };


                for (auto const& sourceParticlesArray : particlesArrays)
                {
                    for (auto const& particle : *sourceParticlesArray)
                    {
                        std::vector<core::Particle<dim>> refinedParticles;
                        auto particleRefinedPos{particle};

                        for (auto iDim = 0u; iDim < dim; ++iDim)
                        {
                            particleRefinedPos.iCell[iDim]
                                = particle.iCell[iDim] * ratio[iDim]
                                  + static_cast<int>(particle.delta[iDim] * ratio[iDim]);
                            particleRefinedPos.delta[iDim]
                                = particle.delta[iDim] * ratio[iDim]
                                  - static_cast<int>(particle.delta[iDim] * ratio[iDim]);
                        }


                        if (isCandidateForSplit_(particleRefinedPos, destinationBox))
                        {
                            split(particleRefinedPos, refinedParticles);


                            // we need to know in which of interior or levelGhostParticlesXXXX
                            // arrays we must put particles

                            bool constexpr putParticlesInCoarseBoundary
                                = splitType == ParticlesDataSplitType::coarseBoundary
                                  || splitType == ParticlesDataSplitType::coarseBoundaryOld
                                  || splitType == ParticlesDataSplitType::coarseBoundaryNew;



                            if (putParticlesInCoarseBoundary)
                            {
                                if constexpr (splitType == ParticlesDataSplitType::coarseBoundary)
                                {
                                    std::copy_if(
                                        std::begin(refinedParticles), std::end(refinedParticles),
                                        std::back_inserter(destCoarseBoundaryParticles), isInDest);
                                }
                                else if constexpr (splitType
                                                   == ParticlesDataSplitType::coarseBoundaryOld)
                                {
                                    //
                                    std::copy_if(std::begin(refinedParticles),
                                                 std::end(refinedParticles),
                                                 std::back_inserter(destCoarseBoundaryOldParticles),
                                                 isInDest);
                                }
                                else //  splitType is coarseBoundaryNew
                                {
                                    std::copy_if(std::begin(refinedParticles),
                                                 std::end(refinedParticles),
                                                 std::back_inserter(destCoarseBoundaryNewParticles),
                                                 isInDest);
                                }
                            }

                            else
                            {
                                std::copy_if(std::begin(refinedParticles),
                                             std::end(refinedParticles),
                                             std::back_inserter(destDomainParticles), isInDest);
                            }
                        } // end is candidate for split
                    }     // end loop on particles
                }         // end loop on source particle arrays
            }             // loop on destination box
        }

        // constexpr int maxCellDistanceFromSplit() const { return std::ceil((interpOrder + 1) *
        // 0.5); }

        SAMRAI::hier::Box getSplitBox(SAMRAI::hier::Box destinationBox) const
        {
            SAMRAI::hier::Box splitBox{destinationBox};
            SAMRAI::tbox::Dimension dimension{dim};
            auto growingVec = SAMRAI::hier::IntVector::getZero(dimension);

            for (auto iDim = 0u; iDim < dim; ++iDim)
            {
                growingVec[iDim] = SplitT::maxCellDistanceFromSplit();
            }
            splitBox.grow(growingVec);

            return splitBox;
        }

        template<typename Particle>
        bool isCandidateForSplit_(Particle const& particle,
                                  SAMRAI::hier::Box const& toFillBox) const
        {
            auto toSplitBox = getSplitBox(toFillBox);
            return isInBox(toSplitBox, particle);
        }
    };
} // namespace amr

} // namespace PHARE

#endif
