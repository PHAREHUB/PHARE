#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_H
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_H

#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/MemoryUtilities.h>

#include <numeric>

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
        , interiorLocalBox_{box}
    {
        interiorLocalBox_.setLower(interiorLocalBox_.lower() - interiorLocalBox_.lower());
        interiorLocalBox_.setUpper(interiorLocalBox_.upper() - box.lower());
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
            SAMRAI::hier::Box const& sourceBox      = pSource->getGhostBox();
            SAMRAI::hier::Box const& destinationBox = getGhostBox();
            SAMRAI::hier::Box intersectionBox{sourceBox * destinationBox};


            if (!intersectionBox.empty())
            {
                copy_(sourceBox, destinationBox, intersectionBox, *pSource);
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

        std::vector<Particle<dim>> specie;

        if (!pOverlap->isOverlapEmpty())
        {
            size_t numberParticles = 0;
            stream >> numberParticles;
            specie.resize(numberParticles);
            stream.unpack(specie.data(), numberParticles);

            SAMRAI::hier::Transformation const& transformation = pOverlap->getTransformation();
            if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
            {
                SAMRAI::hier::BoxContainer const& boxContainer
                    = pOverlap->getDestinationBoxContainer();
                for (auto const& box : boxContainer)
                {
                    SAMRAI::hier::Box unpackBox{box};
                    SAMRAI::hier::Box intersect{unpackBox * getGhostBox() * box};


                    SAMRAI::hier::Box intersectLocalSource{intersect};


                    auto const& ghostWidth = getGhostCellWidth();


                    for (auto const& particle : specie)
                    {
                        if (isInBox_(intersectLocalSource, particle))
                        {
                            auto shiftedParticle{particle};

                            if (isInBox_(interiorLocalBox_, shiftedParticle))
                            {
                                shiftParticle_(ghostWidth, shiftedParticle);
                                interior.push_back(std::move(shiftedParticle));
                            }
                            else
                            {
                                shiftParticle_(ghostWidth, shiftedParticle);
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

        SAMRAI::hier::IntVector shift{sourceBox.lower() - destinationBox.lower() - oneIndex};
        SAMRAI::hier::Box intersectLocalSource{intersectionBox};

        intersectLocalSource.setLower(intersectionBox.lower() - sourceBox.lower());
        intersectLocalSource.setUpper(intersectionBox.upper() - sourceBox.lower());

        copy_(source, shift, intersectLocalSource);
    }




    void copyWithTransform_(SAMRAI::hier::Box const& sourceBox,
                            SAMRAI::hier::Box const& intersectionBox,
                            SAMRAI::hier::Transformation const& transformation,
                            ParticlesData const& source)
    {
        SAMRAI::hier::IntVector particleIndexShift{transformation.getOffset()};
        SAMRAI::hier::Box intersectLocalSource{intersectionBox};

        transformation.inverseTransform(intersectLocalSource);

        SAMRAI::hier::IntVector const one{
            SAMRAI::hier::IntVector::getOne(SAMRAI::tbox::Dimension{dim})};

        SAMRAI::hier::Index oneIndexShift{one};

        intersectLocalSource.setLower(intersectLocalSource.lower() - sourceBox.lower() - one);
        intersectLocalSource.setUpper(intersectLocalSource.upper() - sourceBox.lower() - one);

        particleIndexShift += sourceBox.lower();
        particleIndexShift += one;

        copy_(source, particleIndexShift, intersectLocalSource);
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
                        shiftParticle_(ghostWidth, shiftedParticle);
                        interior.push_back(std::move(shiftedParticle));
                    }
                    else
                    {
                        shiftParticle_(ghostWidth, shiftedParticle);
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

    void pack_(std::vector<Particle<dim>>& buffer, SAMRAI::hier::Box const& intersectionBox,
               SAMRAI::hier::Box const& sourceBox,
               SAMRAI::hier::Transformation const& transformation) const
    {
        SAMRAI::hier::Box localSource{intersectionBox};


        transformation.inverseTransform(localSource);

        SAMRAI::hier::IntVector shift{transformation.getOffset()};


        SAMRAI::hier::IntVector const one{
            SAMRAI::hier::IntVector::getOne(SAMRAI::tbox::Dimension{dim})};

        SAMRAI::hier::Index const oneIndex{
            SAMRAI::hier::IntVector::getOne(SAMRAI::tbox::Dimension{dim})};

        localSource.setLower(localSource.lower() - sourceBox.lower() - oneIndex);
        localSource.setUpper(localSource.upper() - sourceBox.lower() - oneIndex);

        shift += sourceBox.lower();
        shift += one;


        std::array<decltype(interior) const*, 2> particlesArray{&interior, &ghost};

        for (auto const& sourceParticlesArray : particlesArray)
        {
            for (auto const& particle : *sourceParticlesArray)
            {
                if (isInBox_(localSource, particle))
                {
                    auto shiftedParticle{particle};
                    shiftParticle_(shift, shiftedParticle);
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
