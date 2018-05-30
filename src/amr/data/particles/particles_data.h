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
        // interiorLocalBox_.setLower(interiorLocalBox_.lower() - box.lower());
        // interiorLocalBox_.setUpper(interiorLocalBox_.upper() - box.lower());
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




    virtual void unpackStream(SAMRAI::tbox::MessageStream& stream,
                              SAMRAI::hier::BoxOverlap const& overlap) final
    {
        SAMRAI::pdat::CellOverlap const* pOverlap
            = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap);
        TBOX_ASSERT(pOverlap != nullptr);

        std::vector<Particle<dim>> particleArray;

        if (!pOverlap->isOverlapEmpty())
        {
            size_t numberParticles = 0;
            stream >> numberParticles;
            particleArray.resize(numberParticles);
            stream.unpack(particleArray.data(), numberParticles);

            SAMRAI::hier::Transformation const& transformation = pOverlap->getTransformation();
            if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
            {
                SAMRAI::hier::BoxContainer const& boxContainer
                    = pOverlap->getDestinationBoxContainer();
                for (auto const& box : boxContainer)
                {
                    // SAMRAI::hier::Box unpackBox{box};
                    auto const intersect = getGhostBox() * box;
                    // SAMRAI::hier::Box intersectLocalSource{intersect};
                    auto intersectLocalSource = AMRToLocal(intersect, getGhostBox());


                    auto particleShift = AMRToLocal(getGhostBox());


                    for (auto const& particle : particleArray)
                    {
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

    ParticleArray<dim> interior;
    ParticleArray<dim> ghost;
    ParticleArray<dim> coarseToFine;



    // Private
private:
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




    /**
     * @brief copy_ prepares the number of cells by which the particle will be shifted
     * Note that the source and destination boxes are GHOST boxes and the ghost width
     * for particle data is one cell (no particle can exit the domain by more than one cell)
     * Thus, the local index of the lower
     *
     *
     * @param sourceBox is the source box in AMR cell-centered index space
     * @param destinationBox is the destination box in AMR cell-centered index space
     * @param intersectionBox is the intersection between source and destination box, in AMR
     * cell-centered index space
     * @param source is the ParticleData patchData from which particles are to be copied into this.
     */
    void copy_(SAMRAI::hier::Box const& sourceBox, SAMRAI::hier::Box const& destinationBox,
               SAMRAI::hier::Box const& intersectionBox, ParticlesData const& source)
    {
        SAMRAI::hier::Index oneIndex{SAMRAI::hier::IntVector::getOne(SAMRAI::tbox::Dimension{dim})};
        SAMRAI::hier::IntVector shift{sourceBox.lower() - destinationBox.lower() /*- oneIndex*/};
        SAMRAI::hier::Box intersectionLocalSource = AMRToLocal(intersectionBox, sourceBox);
        copy_(source, shift, intersectionLocalSource);
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




    void copy_(ParticlesData const& source, SAMRAI::hier::IntVector const& shiftToDestination,
               SAMRAI::hier::Box const& localSource)
    {
        std::array<decltype(source.interior) const*, 2> particlesArray{&source.interior,
                                                                       &source.ghost};
        auto const& ghostWidth = source.getGhostCellWidth();

        for (auto const& sourceParticlesArray : particlesArray)
        {
            for (auto const& particle : *sourceParticlesArray)
            {
                if (isInBox_(localSource, particle))
                {
                    auto shiftedParticle{particle};

                    shiftParticle_(shiftToDestination, shiftedParticle);

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
