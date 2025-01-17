#ifndef PHARE_PARTICLES_DATA_SPLIT_HPP
#define PHARE_PARTICLES_DATA_SPLIT_HPP


#include "core/def/phare_mpi.hpp"

#include "core/def.hpp"
#include "amr/data/particles/particles_data.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "split.hpp"
#include "core/utilities/constants.hpp"
#include "phare_core.hpp"
#include "amr/amr_constants.hpp"

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/pdat/CellOverlap.h>

#include <functional>


namespace PHARE
{
namespace amr
{
    enum class ParticlesDataSplitType {
        coarseBoundary,
        interior,
        coarseBoundaryOld,
        coarseBoundaryNew
    };


    template<std::size_t interp, typename Particle>
    NO_DISCARD Particle toFineGrid(Particle toFine)
    {
        constexpr auto dim   = Particle::dimension;
        constexpr auto ratio = PHARE::amr::refinementRatio;

        for (size_t iDim = 0; iDim < dim; ++iDim)
        {
            auto fineDelta     = toFine.delta[iDim] * ratio;
            int fineDeltaInt   = static_cast<int>(fineDelta);
            toFine.iCell[iDim] = toFine.iCell[iDim] * ratio + fineDeltaInt;
            toFine.delta[iDim] = fineDelta - fineDeltaInt;
        }

        return toFine;
    }


    template<typename ParticleArray, ParticlesDataSplitType splitType, typename Splitter>
    class ParticlesRefineOperator : public SAMRAI::hier::RefineOperator
    {
    public:
        static constexpr auto dim           = Splitter::dimension;
        static constexpr auto interpOrder   = Splitter::interp_order;
        static constexpr auto nbRefinedPart = Splitter::nbRefinedPart;

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
                            SAMRAI::hier::IntVector const& /*ratio*/) const override
        {
            // For the particles we index them as a CellIndex (for the iCell)
            // ie the particles in the iCell live between lower left node of the iCell
            // and upper right of the same iCell
            auto const& destinationFieldOverlap
                = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(fineOverlap);


            // We then need to get our ParticlesData from the patch
            auto destinationParticlesData = std::dynamic_pointer_cast<ParticlesData<ParticleArray>>(
                destination.getPatchData(destinationComponent));

            auto const sourceParticlesData
                = std::dynamic_pointer_cast<ParticlesData<ParticleArray>>(
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
            refine_(*destinationParticlesData, *sourceParticlesData, destinationFieldOverlap);
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
         * splitting of coarse particles onto the destination patch
         */
        void refine_(ParticlesData<ParticleArray>& destParticlesData,
                     ParticlesData<ParticleArray> const& srcParticlesData,
                     SAMRAI::pdat::CellOverlap const& destFieldOverlap) const
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

            Splitter split;

            // The PatchLevelFillPattern had compute boxes that correspond to the expected filling.
            // In case of a coarseBoundary it will most likely give multiple boxes
            // in case of interior, this will be just one box usually
            for (auto const& destinationBox : destBoxes)
            {
                std::array particlesArrays{&srcInteriorParticles, &srcGhostParticles};
                auto splitBox = getSplitBox(destinationBox);

                auto isInDest = [&destinationBox](auto const& particle) //
                { return isInBox(destinationBox, particle); };


                for (auto const& sourceParticlesArray : particlesArrays)
                {
                    for (auto const& particle : *sourceParticlesArray)
                    {
                        std::array<typename ParticleArray::value_type, nbRefinedPart>
                            refinedParticles;
                        auto particleRefinedPos = toFineGrid<interpOrder>(particle);

                        if (isInBox(splitBox, particleRefinedPos))
                        {
                            split(particleRefinedPos, refinedParticles);


                            // we need to know in which of interior or levelGhostParticlesXXXX
                            // arrays we must put particles

                            bool constexpr putParticlesInCoarseBoundary
                                = splitType == ParticlesDataSplitType::coarseBoundary
                                  || splitType == ParticlesDataSplitType::coarseBoundaryOld
                                  || splitType == ParticlesDataSplitType::coarseBoundaryNew;



                            if constexpr (putParticlesInCoarseBoundary)
                            {
                                if constexpr (splitType == ParticlesDataSplitType::coarseBoundary)
                                {
                                    /*std::cout << "copying " << refinedParticles.size()
                                              << " particles into levelGhost\n";*/
                                    std::copy_if(
                                        std::begin(refinedParticles), std::end(refinedParticles),
                                        std::back_inserter(destCoarseBoundaryParticles), isInDest);
                                }
                                else if constexpr (splitType
                                                   == ParticlesDataSplitType::coarseBoundaryOld)
                                {
                                    /*std::cout << "copying " << refinedParticles.size()
                                              << " particles into levelGhostOld\n";*/
                                    std::copy_if(std::begin(refinedParticles),
                                                 std::end(refinedParticles),
                                                 std::back_inserter(destCoarseBoundaryOldParticles),
                                                 isInDest);
                                }
                                else //  splitType is coarseBoundaryNew
                                {
                                    /*std::cout << "copying " << refinedParticles.size()
                                              << " particles into levelGhostNew\n";*/
                                    std::copy_if(std::begin(refinedParticles),
                                                 std::end(refinedParticles),
                                                 std::back_inserter(destCoarseBoundaryNewParticles),
                                                 isInDest);
                                }
                            }

                            else
                            {
                                /*std::cout << "copying " << refinedParticles.size()
                                          << " particles into domain\n";*/
                                std::copy_if(std::begin(refinedParticles),
                                             std::end(refinedParticles),
                                             std::back_inserter(destDomainParticles), isInDest);
                            }
                        } // end is candidate for split
                    } // end loop on particles
                } // end loop on source particle arrays
            } // loop on destination box
        }


        SAMRAI::hier::Box getSplitBox(SAMRAI::hier::Box const& destinationBox) const
        {
            SAMRAI::hier::Box splitBox{destinationBox};
            SAMRAI::tbox::Dimension dimension{dim};
            auto growingVec = SAMRAI::hier::IntVector::getZero(dimension);

            for (auto iDim = 0u; iDim < dim; ++iDim)
            {
                growingVec[iDim] = Splitter::maxCellDistanceFromSplit();
            }
            splitBox.grow(growingVec);

            return splitBox;
        }
    };

} // namespace amr
} // namespace PHARE


namespace PHARE::amr
{
template<typename ParticleArray, typename Splitter>
struct RefinementParams
{
    using InteriorParticleRefineOp
        = ParticlesRefineOperator<ParticleArray, ParticlesDataSplitType::interior, Splitter>;

    using CoarseToFineRefineOpOld
        = ParticlesRefineOperator<ParticleArray, ParticlesDataSplitType::coarseBoundaryOld,
                                  Splitter>;

    using CoarseToFineRefineOpNew
        = ParticlesRefineOperator<ParticleArray, ParticlesDataSplitType::coarseBoundaryNew,
                                  Splitter>;
};

} // namespace PHARE::amr


#endif
