#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP

#include <numeric>
#include <stdexcept>

#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/MemoryUtilities.h>

#include "core/logger.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particles.hpp"
#include "core/data/particles/particle_array.hpp"
#include "amr/resources_manager/amr_utils.hpp"


namespace PHARE::amr
{
template<std::size_t dim>
inline bool isInBox(SAMRAI::hier::Box const& box, std::array<int, dim> const& iCell)
{
    auto const& lower = box.lower();
    auto const& upper = box.upper();

    if (iCell[0] >= lower(0) && iCell[0] <= upper(0))
    {
        if constexpr (dim > 1)
        {
            if (iCell[1] >= lower(1) && iCell[1] <= upper(1))
            {
                if constexpr (dim > 2)
                {
                    if (iCell[2] >= lower(2) && iCell[2] <= upper(2))
                    {
                        return true;
                    }
                }
                else
                {
                    return true;
                }
            }
        }
        else
        {
            return true;
        }
    }
    return false;
}
} // namespace PHARE::amr

namespace PHARE
{
namespace amr
{
    template<typename Particle>
    inline bool isInBox(SAMRAI::hier::Box const& box, Particle const& particle)
    {
        return isInBox(box, particle.iCell_);
    }


    /** @brief ParticlesData is a concrete SAMRAI::hier::PatchData subclass to store Particle data
     *
     * This class encapsulates particle storage known by the module core, and by being derived
     * from PatchData is compatible with the SAMRAI data management system.
     *
     * A ParticlesData encapsulates **three** different particle arrays:
     *
     * - domainParticles : these particles are those for which iCell is within the physical domain
     * of the patch
     *
     * - patchGhostParticles: these particles are located within the ghost layer around the physical
     * domain of the patch. We call the "ghost layer" the layer of ghostCellWidth just outside the
     * physical domain of the patch, on borders that have neighbors patchs of the same level.
     * All the particles in the ghost layer are exact clones of particles located on a neighbor
     * patch of the same level. The ghost particles are getting here when then exit the neighbor
     * patch, and can enter the patch.
     *
     * - levelGhostParticles: these particles are located in a layer just passed the patch
     * boundaries that also are level boundaries. These particles are getting here when there is a
     * particle refinement from a coarser level
     *
     */
    /**
     * @brief The ParticlesData class
     */
    template<typename ParticleArray>
    class ParticlesData : public SAMRAI::hier::PatchData
    {
        static constexpr auto dim = ParticleArray::dimension;

    public:
        ParticlesData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost)
            : SAMRAI::hier::PatchData::PatchData(box, ghost)
            , pack{&domainParticles, &patchGhostParticles, &levelGhostParticles,
                   &levelGhostParticlesOld, &levelGhostParticlesNew}
            , interiorLocalBox_{AMRToLocal(box, this->getGhostBox())}
        {
        }




        ParticlesData()                     = delete;
        ParticlesData(ParticlesData const&) = delete;
        ParticlesData(ParticlesData&&)      = default;



        ParticlesData& operator=(ParticlesData const&) = delete;




        // SAMRAI interface


        /**
         * @brief copy takes a source PatchData and tries to copy its particles
         * in our particle arrays where the source ghostbox and our ghost box overlap
         *
         *
         * this function just hands the source ghost box, our ghost box and their
         * intersection to a private copy_
         *
         * We follow the procedure suggested by SAMRAI PatchDatas. If by anychance
         * the source PatchData was not a ParticleData like we are, we'd give ourselves
         * to its copy2 function, assuming it can copy its source content into us.
         */
        void copy(SAMRAI::hier::PatchData const& source) override
        {
            PHARE_LOG_SCOPE("ParticlesData::copy");

            TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source);

            // throws if fails
            auto& pSource = dynamic_cast<ParticlesData const&>(source);

            SAMRAI::hier::Box const& sourceGhostBox = pSource.getGhostBox();
            SAMRAI::hier::Box const& myGhostBox     = getGhostBox();
            const SAMRAI::hier::Box intersectionBox{sourceGhostBox * myGhostBox};

            if (!intersectionBox.empty())
            {
                copy_(sourceGhostBox, myGhostBox, intersectionBox, pSource);
            }
        }


        /**
         * @brief our copy2 will be called by a PatchData if a ParticleData was
         * given to be copied into another kind of PatchData. Here we chose that
         * copy2 throws unconditiionnally.
         */
        void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination) const override
        {
            throw std::runtime_error("Cannot cast");
        }


        /**
         * @brief copy with an overlap. Does the copy as the other overload but this time
         * the copy must account for the intersection with the boxes within the overlap
         * The copy is done between the source patch data and myself
         */
        void copy(SAMRAI::hier::PatchData const& source,
                  SAMRAI::hier::BoxOverlap const& overlap) override
        {
            PHARE_LOG_SCOPE("ParticlesData::copy with overlap");

            // casts throw on failure
            auto& pSource  = dynamic_cast<ParticlesData const&>(source);
            auto& pOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap);

            SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
            if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
            {
                SAMRAI::hier::BoxContainer const& boxList = pOverlap.getDestinationBoxContainer();
                for (auto const& overlapBox : boxList)
                {
                    SAMRAI::hier::Box sourceGhostBox = pSource.getGhostBox();
                    SAMRAI::hier::Box myGhostBox     = this->getGhostBox();
                    SAMRAI::hier::Box intersectionBox{sourceGhostBox.getDim()};

                    if (isSameBlock(transformation))
                    {
                        if (offsetIsZero(transformation))
                        {
                            intersectionBox = overlapBox * sourceGhostBox * myGhostBox;

                            if (!intersectionBox.empty())
                            {
                                copy_(sourceGhostBox, myGhostBox, intersectionBox, pSource);
                            }
                        }
                        else
                        {
                            SAMRAI::hier::Box shiftedSourceBox{sourceGhostBox};
                            transformation.transform(shiftedSourceBox);
                            intersectionBox = overlapBox * shiftedSourceBox * myGhostBox;


                            if (!intersectionBox.empty())
                            {
                                copyWithTransform_(sourceGhostBox, intersectionBox, transformation,
                                                   pSource);
                            }
                        }
                    }
                    else
                    {
                        std::runtime_error("Error - multiblock hierarchies not handled");
                    }
                } // end loop over boxes
            }     // end no rotate
            else
            {
                throw std::runtime_error("copy with rotate not implemented");
            }
        }


        void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination,
                   [[maybe_unused]] SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            throw std::runtime_error("Cannot cast");
        }



        bool canEstimateStreamSizeFromBox() const override { return false; }




        std::size_t getDataStreamSize(SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};

            return countNumberParticlesIn_(pOverlap) * ParticleArray::size_of_particle();
        }




        /**
         * @brief packStream is the function that takes particles from our particles arrays
         * that lie in the boxes of the given overlap, and pack them to a stream.
         *
         * Streaming particles means that we have to take particles with iCell on a local source
         * index space , communicate them, and load them at destination with iCell on a destination
         * local index space. To do that we need to:
         *
         * 1- translate source iCell to source AMR index space
         * 2- Apply the offset to shift this AMR index on top of the destination cells
         * 3- pack and communicate particles
         * 4- move back iCell from the shifted AMR index space to the local destination index space
         *
         * Note that step 2 could be done upon reception of the pack, we chose to do it before.
         *
         */
        void packStream(SAMRAI::tbox::MessageStream& stream,
                        SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            PHARE_LOG_SCOPE("packStream");

            auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};

            ParticleArray specie;

            if (pOverlap.isOverlapEmpty())
            {
                constexpr std::size_t zero = 0;
                stream << zero;
            }
            else
            {
                SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
                if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
                {
                    SAMRAI::hier::BoxContainer const& boxContainer
                        = pOverlap.getDestinationBoxContainer();

                    auto const& sourceGhostBox = getGhostBox();

                    // sourceBox + offset = source on destination
                    // we are given boxes in the Overlap in destination
                    // index space. And we want to select all particles
                    // in the ghost source box that lie in this overlapBox
                    // we thus need to first shift the sourceGhostBox to the
                    // destination index space so that its cells (partly) overlap the one
                    // of the given overlap boxes.
                    // Then pack_ will take all particles which iCell, shifted by the
                    // transformation offset onto the overlapBox index space,
                    // lie in the overlap box.
                    SAMRAI::hier::Box transformedSource{sourceGhostBox};
                    transformation.transform(transformedSource);

                    for (auto const& overlapBox : boxContainer)
                    {
                        SAMRAI::hier::Box intersectionBox{transformedSource * overlapBox};

                        pack_(specie, intersectionBox, sourceGhostBox, transformation);
                    }
                }
                else
                {
                    throw std::runtime_error("Error - rotations not handled in PHARE");
                }
                stream << specie.size();
                stream.growBufferAsNeeded();
                this->pack_(stream, specie);
            }
        }


        void pack_(SAMRAI::tbox::MessageStream& stream, ParticleArray const& specie) const
        {
            if constexpr (ParticleArray::is_contiguous)
                std::apply(
                    [&](auto const&... container) {
                        ((stream.pack(container.data(), specie.size())), ...);
                    },
                    specie.as_tuple());
            else
                stream.pack(specie.data(), specie.size());
        }


        /**
         * @brief unpackStream is the function that unpacks a stream of particles to our particle
         * arrays.
         *
         * We get a stream and an overlap. The overlap contains boxes where to put particles and
         * transformation from source to destination AMR indexes.
         *
         * By convention chosen in patckStream, packed particles have their iCell in our AMR index
         * space. This means that before putting them into our local arrays, we need to apply
         * AMRToLocal() to get the proper shift to apply to them
         *
         */
        void unpackStream(SAMRAI::tbox::MessageStream& stream,
                          SAMRAI::hier::BoxOverlap const& overlap) override
        {
            PHARE_LOG_SCOPE("unpackStream");

            auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};

            if (!pOverlap.isOverlapEmpty())
            {
                // unpack particles into a particle array
                std::size_t numberParticles = 0;
                stream >> numberParticles;
                ParticleArray particleArray(numberParticles);
                unpack_(stream, particleArray);

                // ok now our goal is to put the particles we have just unpacked
                // into the particleData and in the proper particleArray : interior or ghost

                SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
                if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
                {
                    // we loop over all boxes in the overlap
                    // we have to first take the intersection of each of these boxes
                    // with our ghostBox. This is where unpacked particles should go.

                    SAMRAI::hier::BoxContainer const& overlapBoxes
                        = pOverlap.getDestinationBoxContainer();

                    auto myBox      = getBox();
                    auto myGhostBox = getGhostBox();

                    for (auto const& overlapBox : overlapBoxes)
                    {
                        // our goal here is :
                        // 1/ to check if each particle is in the intersect of the overlap boxes
                        // and our ghostBox 2/ if yes, check if these particles should go within the
                        // interior array or ghost array
                        auto const intersect = getGhostBox() * overlapBox;

                        for (std::size_t p_idx = 0; p_idx < particleArray.size(); ++p_idx)
                        {
                            if (isInBox(intersect, particleArray.iCell(p_idx)))
                            {
                                if (isInBox(myBox, particleArray.iCell(p_idx)))
                                {
                                    domainParticles.push_back(particleArray[p_idx]);
                                }
                                else
                                {
                                    patchGhostParticles.push_back(particleArray[p_idx]);
                                }
                            }
                        } // end species loop
                    }     // end box loop
                }         // end no rotation
            }             // end overlap not empty
        }

        void unpack_(SAMRAI::tbox::MessageStream& stream, ParticleArray& specie) const
        {
            if constexpr (ParticleArray::is_contiguous)
                std::apply(
                    [&](auto&... container) {
                        ((stream.unpack(container.data(), specie.size())), ...);
                    },
                    specie.as_tuple());
            else
                stream.unpack(specie.data(), specie.size());
        }


        core::ParticlesPack<ParticleArray>* getPointer() { return &pack; }



        // Core interface
        // these particles arrays are public because core module is free to use
        // them easily
        ParticleArray domainParticles;
        ParticleArray patchGhostParticles;

        ParticleArray levelGhostParticles;

        ParticleArray levelGhostParticlesOld;
        ParticleArray levelGhostParticlesNew;

        core::ParticlesPack<ParticleArray> pack;



    private:
        //! interiorLocalBox_ is the box, in local index space, that goes from the first to the last
        //! cell in our patch physical domain, i.e. "from dual physical start index to dual physical
        //! end index"
        SAMRAI::hier::Box interiorLocalBox_;



        void copy_(SAMRAI::hier::Box const& /*sourceGhostBox*/,
                   SAMRAI::hier::Box const& /*destinationGhostBox*/,
                   SAMRAI::hier::Box const& intersectionBox, ParticlesData const& sourceData)
        {
            std::array<ParticleArray const*, 2> particlesArrays{&sourceData.domainParticles,
                                                                &sourceData.patchGhostParticles};

            auto myDomainBox = this->getBox();

            // for each particles in the source ghost and domain particle arrays
            // we check if it is in the intersectionBox
            // if it is, is it in my domain box ?
            //      - if so, let's add it to my domain particle array
            //      - if not, let's add it to my ghost particle array
            for (auto const& sourceParticlesArray : particlesArrays)
            {
                auto& particleArray = *sourceParticlesArray;
                for (std::size_t p_idx = 0; p_idx < particleArray.size(); ++p_idx)
                {
                    if (isInBox(intersectionBox, particleArray.iCell(p_idx)))
                    {
                        if (isInBox(myDomainBox, particleArray.iCell(p_idx)))
                        {
                            domainParticles.push_back(particleArray[p_idx]);
                        }
                        else
                        {
                            patchGhostParticles.push_back(particleArray[p_idx]);
                        }
                    }
                }
            }
        }




        void copyWithTransform_(SAMRAI::hier::Box const& /*sourceGhostBox*/,
                                SAMRAI::hier::Box const& intersectionBox,
                                SAMRAI::hier::Transformation const& transformation,
                                ParticlesData const& sourceData)
        {
            std::array<ParticleArray const*, 2> particlesArrays{&sourceData.domainParticles,
                                                                &sourceData.patchGhostParticles};

            auto myDomainBox = this->getBox();

            auto offset = transformation.getOffset();

            for (auto const& sourceParticlesArray : particlesArrays)
            {
                auto& particleArray = *sourceParticlesArray;
                for (std::size_t p_idx = 0; p_idx < particleArray.size(); ++p_idx)
                {
                    // the particle is only copied if it is in the intersectionBox
                    // but before its iCell must be shifted by the transformation offset

                    auto newParticleiCell{particleArray.iCell(p_idx)};
                    for (auto iDir = 0u; iDir < newParticleiCell.size(); ++iDir)
                        newParticleiCell[iDir] += offset[iDir];

                    if (isInBox(intersectionBox, newParticleiCell))
                    {
                        // now we now the particle is in the intersection
                        // we need to know whether it is in the domain part of that
                        // intersection. If it is not, then it must be in the ghost part

                        auto add_update = [&](auto& array) {
                            array.push_back(particleArray[p_idx]);
                            array.iCell(array.size() - 1) = newParticleiCell;
                        };

                        if (isInBox(myDomainBox, newParticleiCell))
                            add_update(domainParticles);
                        else
                            add_update(patchGhostParticles);
                    }
                }
            }


            // SAMRAI::hier::Box localSourceSelectionBox = AMRToLocal(intersectionBox,
            // sourceGhostBox);

            // we shift it back the box on top of source AMR indexes
            // transformation.inverseTransform(localSourceSelectionBox);

            // copy_(sourceData, particleShift, localSourceSelectionBox);
        }




        /**
         * @brief countNumberParticlesIn_ counts the number of particles that lie
         * within the boxes of an overlap. This function count both patchGhost and
         * domain particles since both could be streamed and we want an upperbound
         * on the number of bytes that could be streamed.
         */
        std::size_t countNumberParticlesIn_(SAMRAI::pdat::CellOverlap const& overlap) const
        {
            std::size_t numberParticles = 0;

            if (overlap.isOverlapEmpty())
            {
                return numberParticles;
            }

            auto const& overlapBoxes = overlap.getDestinationBoxContainer();

            for (auto const& overlapBox : overlapBoxes)
            {
                // we are given boxes from the overlap
                // we want to know how many of our local particles
                // lie in that overlap. Overlap is given in the destination
                // index space (see overlap documentation)
                // so we need to transform that overlap box into our box index space.
                // Since source index space + offset = destination indexspace
                // we need to apply an inverseTransform to the overlapBox.
                // then we intersect it with our Box and count how many of domain particles
                // our inside that intersection.
                SAMRAI::hier::Box shiftedOverlapBox{overlapBox};
                SAMRAI::hier::Transformation const& transformation = overlap.getTransformation();
                transformation.inverseTransform(shiftedOverlapBox);
                SAMRAI::hier::Box intersectionBox{shiftedOverlapBox * getGhostBox()};

                numberParticles += countNumberParticlesIn_(intersectionBox);
            }
            return numberParticles;
        }


        /**
         * @brief countNumberParticlesIn_ returns the number of interior particles within a given
         * box
         *
         * the box given is in AMR index space so the function first needs to put it in
         * local indexing relative to the domain box
         */
        std::size_t countNumberParticlesIn_(SAMRAI::hier::Box const& box) const
        {
            std::size_t numberParticles{0};

            for (std::size_t i = 0; i < domainParticles.size(); ++i)
                if (isInBox(box, domainParticles.iCell(i)))
                    ++numberParticles;

            return numberParticles;
        }




        void pack_(ParticleArray& buffer, SAMRAI::hier::Box const& intersectionBox,
                   [[maybe_unused]] SAMRAI::hier::Box const& sourceBox,
                   SAMRAI::hier::Transformation const& transformation) const
        {
            std::array<ParticleArray const*, 2> particlesArrays{&domainParticles,
                                                                &patchGhostParticles};

            for (auto const& sourceParticlesArray : particlesArrays)
            {
                auto& particleArray = *sourceParticlesArray;
                for (std::size_t p_idx = 0; p_idx < particleArray.size(); ++p_idx)
                {
                    auto shiftediCell{particleArray.iCell(p_idx)};
                    auto offset = transformation.getOffset();

                    for (auto i = 0u; i < dim; ++i)
                        shiftediCell[i] += offset[i];

                    if (isInBox(intersectionBox, shiftediCell))
                    {
                        buffer.push_back(particleArray[p_idx]);
                        buffer.iCell(buffer.size() - 1) = shiftediCell;
                    }
                }
            }
        }
    };

} // namespace amr
} // namespace PHARE

#endif
