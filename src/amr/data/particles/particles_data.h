#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_H
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_H

#include <numeric>
#include <stdexcept>

#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/MemoryUtilities.h>


#include "data/particles/particle.h"
#include "data/particles/particle_array.h"

namespace PHARE
{
template<std::size_t dim>
class ParticlesData : public SAMRAI::hier::PatchData
{
public:
    ParticlesData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost)
        : SAMRAI::hier::PatchData::PatchData(box, ghost)
        , interiorLocalBox_{SAMRAI::hier::Index{box.getDim(), 0} + ghost, box.upper() - box.lower(),
                            box.getBlockId()}
    {
    }




    ParticlesData()                     = delete;
    ParticlesData(ParticlesData const&) = delete;
    ParticlesData(ParticlesData&&)      = default;



    ParticlesData& operator=(ParticlesData const&) = delete;




    // SAMRAI interface

    virtual void copy(SAMRAI::hier::PatchData const& source) final
    {
        TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source);

        const ParticlesData* pSource = dynamic_cast<const ParticlesData*>(&source);
        if (pSource != nullptr)
        {
            SAMRAI::hier::Box const& sourceGhostBox = pSource->getGhostBox();
            SAMRAI::hier::Box const& myGhostBox     = getGhostBox();
            SAMRAI::hier::Box intersectionBox{sourceGhostBox * myGhostBox};

            if (!intersectionBox.empty())
            {
                copy_(sourceGhostBox, myGhostBox, intersectionBox, *pSource);
            }
        }
        else
        {
            source.copy2(*this);
        }
    }




    virtual void copy2(SAMRAI::hier::PatchData& destination) const final
    {
        throw std::runtime_error("Cannot cast");
    }




    bool offsetIsZero(SAMRAI::hier::Transformation const& transformation)
    {
        auto const& offset = transformation.getOffset();
        auto dimension     = offset.getDim();
        return transformation.getOffset() == SAMRAI::hier::IntVector::getZero(dimension);
    }




    bool isSameBlock(SAMRAI::hier::Transformation const& transformation)
    {
        return transformation.getBeginBlock() == transformation.getEndBlock();
    }




    virtual void copy(SAMRAI::hier::PatchData const& source,
                      SAMRAI::hier::BoxOverlap const& overlap) final
    {
        const ParticlesData* pSource = dynamic_cast<const ParticlesData*>(&source);
        const SAMRAI::pdat::CellOverlap* pOverlap
            = dynamic_cast<const SAMRAI::pdat::CellOverlap*>(&overlap);


        if ((pSource != nullptr) && (pOverlap != nullptr))
        {
            SAMRAI::hier::Transformation const& transformation = pOverlap->getTransformation();
            if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
            {
                SAMRAI::hier::BoxContainer const& boxList = pOverlap->getDestinationBoxContainer();
                for (auto const& overlapBox : boxList)
                {
                    SAMRAI::hier::Box sourceBox      = pSource->getGhostBox();
                    SAMRAI::hier::Box destinationBox = this->getGhostBox();
                    SAMRAI::hier::Box intersectionBox{sourceBox.getDim()};

                    if (isSameBlock(transformation))
                    {
                        if (offsetIsZero(transformation))
                        {
                            intersectionBox = overlapBox * sourceBox * destinationBox;

                            if (!intersectionBox.empty())
                            {
                                copy_(sourceBox, destinationBox, intersectionBox, *pSource);
                            }
                        }
                        else
                        {
                            SAMRAI::hier::Box shiftedSourceBox{sourceBox};
                            transformation.transform(shiftedSourceBox);
                            intersectionBox = overlapBox * shiftedSourceBox * destinationBox;

                            if (!intersectionBox.empty())
                            {
                                copyWithTransform_(sourceBox, intersectionBox, transformation,
                                                   *pSource);
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
        else
        {
            source.copy2(*this, overlap);
        }
    }




    virtual void copy2(SAMRAI::hier::PatchData& destination,
                       SAMRAI::hier::BoxOverlap const& overlap) const final
    {
        throw std::runtime_error("Cannot cast");
    }




    virtual bool canEstimateStreamSizeFromBox() const final { return false; }




    virtual size_t getDataStreamSize(SAMRAI::hier::BoxOverlap const& overlap) const final
    {
        SAMRAI::pdat::CellOverlap const* pOverlap{
            dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap)};

        std::size_t numberParticles = countNumberParticlesIn_(*pOverlap);

        return SAMRAI::tbox::MemoryUtilities::align(numberParticles * sizeof(ParticleArray<dim>));
    }




    /**
     * @brief packStream is the function that takes particles from our particles arrays
     * that lie in the boxes of the given overlap, and pack them to a stream.
     *
     * Streaming particles means that we have to take particles with iCell on a local source index
     * space , communicate them, and load them at destination with iCell on a destination local
     * index space. To do that we need to:
     *
     * 1- translate source iCell to source AMR index space
     * 2- Apply the offset to shift this AMR index on top of the destination cells
     * 3- pack and communicate particles
     * 4- move back iCell from the shifted AMR index space to the local destination index space
     *
     * Note that step 2 could be done upon reception of the pack, we chose to do it before.
     *
     * Thus, by convention, the local iCell of the packed particles is translated to the
     * destination AMR index space, i.e. it is moved to AMR space and shifted by the
     * offset transformation given in the overlap.
     *
     * example : say we have an AMR domain [0,15] with two patches P1[0,5] and P2[10,15]
     * with periodic boundaries. AMR index 15 is thus equivalent to AMR index -1, and 16 to 0
     *
     * Say we have one ghost cell.
     *
     * A particle on patch P2 with local iCell == 7 is thus on AMR index 16
     * we want to stream it to P1.
     *
     * moving local iCell 7 to AMR means iCell = 16
     * applying offset P2-->P1 means iCell becomes 0 (16 - offset_P2P1=16)
     *
     * at this point we pack it, SAMRAI communicates it to destination
     * in unpack we get iCell = 0
     *
     * then we move it from AMR to P1 local means iCell becomes 1
     * at this point with iCell=1 we know the particle should be placed into the interior particle
     * buffer
     *
     */
    virtual void packStream(SAMRAI::tbox::MessageStream& stream,
                            SAMRAI::hier::BoxOverlap const& overlap) const final
    {
        SAMRAI::pdat::CellOverlap const* pOverlap{
            dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap)};

        TBOX_ASSERT(pOverlap != nullptr);

        std::vector<Particle<dim>> specie;

        if (pOverlap->isOverlapEmpty())
        {
            constexpr std::size_t zero = 0;
            stream << zero;
        }
        else
        {
            SAMRAI::hier::Transformation const& transformation = pOverlap->getTransformation();
            if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
            {
                SAMRAI::hier::BoxContainer const& boxContainer
                    = pOverlap->getDestinationBoxContainer();
                for (auto const& destinationBox : boxContainer)
                {
                    auto const& sourceBox = getGhostBox();

                    SAMRAI::hier::Box transformedSource{sourceBox};
                    transformation.transform(transformedSource);

                    SAMRAI::hier::Box intersectionBox{transformedSource * destinationBox};

                    pack_(specie, intersectionBox, sourceBox, transformation);
                }
            }
            else
            {
                throw std::runtime_error("Error - rotations not handled in PHARE");
            }
            stream << specie.size();
            stream.pack(specie.data(), specie.size());
        }
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
    virtual void unpackStream(SAMRAI::tbox::MessageStream& stream,
                              SAMRAI::hier::BoxOverlap const& overlap) final
    {
        SAMRAI::pdat::CellOverlap const* pOverlap
            = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap);
        TBOX_ASSERT(pOverlap != nullptr);

        std::vector<Particle<dim>> particleArray;

        if (!pOverlap->isOverlapEmpty())
        {
            // unpack particles into a particle array
            size_t numberParticles = 0;
            stream >> numberParticles;
            particleArray.resize(numberParticles);
            stream.unpack(particleArray.data(), numberParticles);

            // ok now our goal is to put the particles we have just unpacked
            // into the particleData and in the proper particleArray : interior or ghost

            SAMRAI::hier::Transformation const& transformation = pOverlap->getTransformation();
            if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
            {
                // we loop over all boxes in the overlap
                // we have to first take the intersection of each of these boxes
                // with our ghostBox. This is where unpacked particles should go.

                SAMRAI::hier::BoxContainer const& destinationBoxes
                    = pOverlap->getDestinationBoxContainer();


                for (auto const& box : destinationBoxes)
                {
                    // our goal here is :
                    // 1/ to check if the each particle is in the intersect of the overlap boxes and
                    // our ghostBox 2/ if yes, check if these particles should go within the
                    // interior array or ghost array

                    // first we need to calculate the intersect between our ghostBox and the
                    // overlapBox these and the resulting intersect are in AMR index space unpacked
                    // particles have their iCell in our AMR index space. we could :
                    //  1- compare their iCell with the intersect right away, if within shit them to
                    //  local indexing and push them in proper array
                    // or
                    // 2- translate the intersect from  AMR to local indexing, shift the particles
                    // from AMR to local, and do the comparison in local index space and then push
                    // them in the proper array

                    auto const intersect = getGhostBox() * box;


                    // we chose first option, so get the intersect in local index space
                    // and get the vector to shift a particle from our AMR space to local index
                    // space
                    auto intersectLocalSource = AMRToLocal(intersect, getGhostBox());
                    auto particleShift        = AMRToLocal(getGhostBox());

                    for (auto const& particle : particleArray)
                    {
                        // shift the particle to local index space
                        // and if it is in intersection, decide in which array to push it.
                        auto shiftedParticle{particle};
                        shiftParticle_(particleShift, shiftedParticle);

                        if (isInBox_(intersectLocalSource, shiftedParticle))
                        {
                            if (isInBox_(interiorLocalBox_, shiftedParticle))
                            {
                                interior.push_back(std::move(shiftedParticle));
                            }
                            else
                            {
                                ghost.push_back(std::move(shiftedParticle));
                            }
                        }
                    } // end species loop
                }     // end box loop
            }         // end no rotation
        }             // end overlap not empty
    }


    // Core interface
    // these particles arrays are public because core module is free to use
    // them easily
    ParticleArray<dim> interior;
    ParticleArray<dim> ghost;
    ParticleArray<dim> coarseToFine;



    // Private
private:
    //! interiorLocalBox_ is the box, in local index space, that goes from the first to the last
    //! cell in our patch physical domain, i.e. "from dual physical start index to dual physical end
    //! index"
    SAMRAI::hier::Box interiorLocalBox_;


    void AMRToLocal(SAMRAI::hier::Box& AMRBox, SAMRAI::hier::Box const& referenceAMR)
    {
        AMRBox.setLower(AMRBox.lower() - referenceAMR.lower());
        AMRBox.setUpper(AMRBox.upper() - referenceAMR.lower());
    }




    SAMRAI::hier::Box AMRToLocal(SAMRAI::hier::Box const& AMRBox,
                                 SAMRAI::hier::Box const& referenceAMR) const
    {
        SAMRAI::hier::Box localBox{AMRBox};
        localBox.setLower(AMRBox.lower() - referenceAMR.lower());
        localBox.setUpper(AMRBox.upper() - referenceAMR.lower());
        return localBox;
    }




    void copy_(SAMRAI::hier::Box const& sourceGhostBox,
               SAMRAI::hier::Box const& destinationGhostBox,
               SAMRAI::hier::Box const& intersectionBox, ParticlesData const& sourceData)
    {
        // at this point indersectionBox is in AMR index space
        // we need to translate it to local source index space to copy_
        SAMRAI::hier::Box intersectionLocalSource = AMRToLocal(intersectionBox, sourceGhostBox);

        // copy_ also needs the vector to shift particles from source local indexing to destination
        // local indexing
        auto particleShift = particleCellShift(sourceGhostBox, destinationGhostBox);

        copy_(sourceData, particleShift, intersectionLocalSource);
    }




    SAMRAI::hier::IntVector localToAMR(SAMRAI::hier::Box const& referenceAMRBox)
    {
        return SAMRAI::hier::IntVector{referenceAMRBox.lower()};
    }




    SAMRAI::hier::IntVector AMRToLocal(SAMRAI::hier::Box const& referenceAMRBox)
    {
        SAMRAI::hier::Index zero{SAMRAI::tbox::Dimension{dim}, 0};
        return SAMRAI::hier::IntVector{zero - referenceAMRBox.lower()};
    }


    /**
     * @brief particleCellShift returns the vector to shift particle iCell from local source
     * indexing to destination source indexing when there is no offset in AMR index between source
     * and destination boxes
     */
    SAMRAI::hier::IntVector particleCellShift(SAMRAI::hier::Box const& sourceGhostBox,
                                              SAMRAI::hier::Box const& destGhostBox)
    {
        auto particleCellShift = localToAMR(sourceGhostBox);
        particleCellShift += AMRToLocal(this->getGhostBox());
        return particleCellShift;
    }


    /**
     * @brief particleCellShift returns the vector to shift particle iCell from local source
     * indexing to destination source indexing with an offset in AMR indexing between source and
     * destination boxes
     */
    SAMRAI::hier::IntVector particleCellShift(SAMRAI::hier::Box const& sourceGhostBox,
                                              SAMRAI::hier::Transformation const& transformation,
                                              SAMRAI::hier::Box const& destGhostBox)
    {
        auto particleCellShift = localToAMR(sourceGhostBox);
        particleCellShift += transformation.getOffset();
        particleCellShift += AMRToLocal(this->getGhostBox());
        return particleCellShift;
    }




    void copyWithTransform_(SAMRAI::hier::Box const& sourceGhostBox,
                            SAMRAI::hier::Box const& intersectionBox,
                            SAMRAI::hier::Transformation const& transformation,
                            ParticlesData const& sourceData)
    {
        auto particleShift = particleCellShift(sourceGhostBox, transformation, this->getGhostBox());

        SAMRAI::hier::Box localSourceSelectionBox = AMRToLocal(intersectionBox, sourceGhostBox);

        // we shift it back the box on top of source AMR indexes
        transformation.inverseTransform(localSourceSelectionBox);

        copy_(sourceData, particleShift, localSourceSelectionBox);
    }



    /**
     * @brief copy_ this overload takes the source patch data, the box in local source indexing
     * where to select particles to copy, and the vector to shift their iCell to our (destination)
     * local indexing and copy them to our local particle buffers.
     */
    void copy_(ParticlesData const& sourceData,
               SAMRAI::hier::IntVector const& shiftParticleCellToDest,
               SAMRAI::hier::Box const& localSourceSelectionBox)
    {
        std::array<decltype(sourceData.interior) const*, 2> particlesArrays{&sourceData.interior,
                                                                            &sourceData.ghost};

        // loop over both interior and ghost source particle buffers
        // for each of their particles, check they should be copied (they should be intersection)
        // and if they do, check into which of our particle buffers it should be copied into.
        // beware the before checking wether it is inside our domain or not, we must
        // translate their iCell from local source indexing to our local indexing.
        for (auto const& sourceParticlesArray : particlesArrays)
        {
            for (auto const& particle : *sourceParticlesArray)
            {
                if (isInBox_(localSourceSelectionBox, particle))
                {
                    auto shiftedParticle{particle};

                    shiftParticle_(shiftParticleCellToDest, shiftedParticle);

                    // TODO isInBox existe déja on peut pas l'utiliser ?
                    if (isInBox_(interiorLocalBox_, shiftedParticle))
                    {
                        // shiftParticle_(ghostWidth, shiftedParticle);
                        interior.push_back(std::move(shiftedParticle));
                    }
                    else
                    {
                        // shiftParticle_(ghostWidth, shiftedParticle);
                        ghost.push_back(std::move(shiftedParticle));
                    }
                }
            }
        }
    }




    bool isInBox_(SAMRAI::hier::Box const& box, Particle<dim> const& particle) const;




    void shiftParticle_(SAMRAI::hier::IntVector const& shift, Particle<dim>& particle) const;




    std::size_t countNumberParticlesIn_(SAMRAI::pdat::CellOverlap const& overlap) const
    {
        std::size_t numberParticles = 0;

        if (overlap.isOverlapEmpty())
        {
            return numberParticles;
        }

        auto const& boxes = overlap.getDestinationBoxContainer();

        for (auto const& box : boxes)
        {
            SAMRAI::hier::Box finalBox{box};
            SAMRAI::hier::Transformation const& transformation = overlap.getTransformation();
            transformation.transform(finalBox);
            SAMRAI::hier::Box intersectionBox{finalBox * getBox()};

            numberParticles += countNumberParticlesIn_(intersectionBox);
        }


        return numberParticles;
    }

    std::size_t countNumberParticlesIn_(SAMRAI::hier::Box const& box) const
    {
        std::size_t numberParticles{0};

        SAMRAI::hier::Box localSource{box};

        auto const& sourceBox = getBox();

        localSource.setLower(box.lower() - sourceBox.lower());
        localSource.setUpper(box.upper() - sourceBox.upper());

        for (auto const& particle : interior)
        {
            if (isInBox_(localSource, particle))
            {
                ++numberParticles;
            }
        }

        return numberParticles;
    }


    /**
     * @brief localToShiftedAMR returns the vector to shift a particle cell from a local
     * index to an AMR offseted index. This function is to be used to shift particles cells
     * before packing them to the stream.
     * At the unpack stage, when the destination box is known, one will have to
     * call AMRToLocal() to finish the index shift to get particles cells into destination local
     * index space
     */
    SAMRAI::hier::IntVector
    localToShiftedAMR(SAMRAI::hier::Box const& sourceGhostBox,
                      SAMRAI::hier::Transformation const& transformation) const
    {
        SAMRAI::hier::IntVector shift{sourceGhostBox.lower()}; // move to AMR space
        shift += transformation.getOffset(); // translate AMR indexes with the transfo offset
        return shift;
    }




    void pack_(std::vector<Particle<dim>>& buffer, SAMRAI::hier::Box const& intersectionBox,
               SAMRAI::hier::Box const& sourceBox,
               SAMRAI::hier::Transformation const& transformation) const
    {
        auto particleShift = localToShiftedAMR(sourceBox, transformation);


        // take the intersectionBox, which is in AMR space aligned with the destination
        // shift its AMR index to the source location
        // and move the shifted AMR index to local index space
        auto localSelectionSourceBox{intersectionBox};
        transformation.inverseTransform(localSelectionSourceBox);
        localSelectionSourceBox = AMRToLocal(localSelectionSourceBox, sourceBox);

        // now it is possible to compare localSelectionSourceBox indexing and
        // particle.iCell on source patch data to select those in the intersection.

        // only the interior and ghost particles are to be copied
        // coarseToFine particles should not be copied.
        std::array<decltype(interior) const*, 2> particlesArrays{&interior, &ghost};

        for (auto const& sourceParticlesArray : particlesArrays)
        {
            for (auto const& particle : *sourceParticlesArray)
            {
                if (isInBox_(localSelectionSourceBox, particle))
                {
                    auto shiftedParticle{particle};
                    shiftParticle_(particleShift, shiftedParticle);
                    buffer.push_back(shiftedParticle);
                }
            }
        }
    }
};




template<>
bool ParticlesData<1>::isInBox_(SAMRAI::hier::Box const& box, Particle<1> const& particle) const
{
    auto const& iCell = particle.iCell;

    auto const& lower = box.lower();
    auto const& upper = box.upper();

    if (iCell[0] >= lower(0) && iCell[0] <= upper(0))
    {
        return true;
    }
    return false;
}

template<>
bool ParticlesData<2>::isInBox_(SAMRAI::hier::Box const& box, Particle<2> const& particle) const
{
    auto const& iCell = particle.iCell;

    auto const& lower = box.lower();
    auto const& upper = box.upper();


    if (iCell[0] >= lower(0) && iCell[0] <= upper(0))
    {
        if (iCell[1] >= lower(1) && iCell[1] <= upper(1))
        {
            return true;
        }
    }
    return false;
}

template<>
bool ParticlesData<3>::isInBox_(SAMRAI::hier::Box const& box, Particle<3> const& particle) const
{
    auto const& iCell = particle.iCell;

    auto const& lower = box.lower();
    auto const& upper = box.upper();


    if (iCell[0] >= lower(0) && iCell[0] <= upper(0))
    {
        if (iCell[1] >= lower(1) && iCell[1] <= upper(1))
        {
            if (iCell[2] >= lower(2) && iCell[2] <= upper(2))
            {
                return true;
            }
        }
    }
    return false;
}

template<>
void ParticlesData<1>::shiftParticle_(SAMRAI::hier::IntVector const& shift,
                                      Particle<1>& particle) const
{
    particle.iCell[0] += shift[0];
}

template<>
void ParticlesData<2>::shiftParticle_(SAMRAI::hier::IntVector const& shift,
                                      Particle<2>& particle) const
{
    particle.iCell[0] += shift[0];
    particle.iCell[1] += shift[1];
}
template<>
void ParticlesData<3>::shiftParticle_(SAMRAI::hier::IntVector const& shift,
                                      Particle<3>& particle) const
{
    particle.iCell[0] += shift[0];
    particle.iCell[1] += shift[1];
    particle.iCell[2] += shift[2];
}

} // namespace PHARE

#endif
