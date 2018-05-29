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
        , interiorBox_{box}
    {
        interiorBox_.setLower(interiorBox_.lower() - interiorBox_.lower());
        interiorBox_.setUpper(interiorBox_.upper() - box.lower());
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
            SAMRAI::hier::Box const& sourceBox      = pSource->getBox();
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
                for (auto const& box : boxList)
                {
                    SAMRAI::hier::Box sourceBox      = pSource->getBox();
                    SAMRAI::hier::Box destinationBox = this->getGhostBox();
                    SAMRAI::hier::Box intersectionBox{sourceBox.getDim()};

                    if (transformation.getOffset() == SAMRAI::hier::IntVector::getZero(box.getDim())
                        && transformation.getBeginBlock() == transformation.getEndBlock())
                    {
                        intersectionBox = box * sourceBox * destinationBox;

                        if (!intersectionBox.empty())
                        {
                            copy_(sourceBox, destinationBox, intersectionBox, *pSource);
                        }
                    }
                    else
                    {
                        SAMRAI::hier::Box transformedSource{sourceBox};
                        transformation.transform(transformedSource);
                        intersectionBox = box * transformedSource * destinationBox;

                        if (!intersectionBox.empty())
                        {
                            copyWithTransform_(sourceBox, transformedSource, destinationBox,
                                               intersectionBox, transformation, *pSource);
                        }
                    }
                }
            }
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
                    auto const& sourceBox = getBox();

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

                    auto const& shift = transformation.getOffset();


                    SAMRAI::hier::Box intersectLocalSource{intersect};

                    transformation.inverseTransform(intersectLocalSource);

                    for (auto const& particle : specie)
                    {
                        if (isInBox_(intersectLocalSource, particle))
                        {
                            auto shiftedParticle{particle};
                            shiftParticle_(shift, shiftedParticle);

                            if (isInBox_(interiorBox_, shiftedParticle))
                            {
                                interior.push_back(std::move(shiftedParticle));
                            }
                            else
                            {
                                ghost.push_back(std::move(shiftedParticle));
                            }
                        }
                    }
                }
            }
        }
    }


    // Core interface

    ParticleArray<dim> interior;
    ParticleArray<dim> ghost;
    ParticleArray<dim> incoming;



    // Private
private:
    SAMRAI::hier::Box interiorBox_;

    void copy_(SAMRAI::hier::Box const& sourceBox, SAMRAI::hier::Box const& destinationBox,
               SAMRAI::hier::Box const& intersectionBox, ParticlesData const& source)
    {
        SAMRAI::hier::IntVector shift{sourceBox.lower() - destinationBox.lower()
                                      - getGhostCellWidth()};
        SAMRAI::hier::Box intersectLocalSource{intersectionBox};

        intersectLocalSource.setLower(intersectionBox.lower() - sourceBox.lower());
        intersectLocalSource.setUpper(intersectionBox.upper() - sourceBox.lower());

        copy_(source, shift, intersectLocalSource);
    }

    void copyWithTransform_(SAMRAI::hier::Box const& sourceBox,
                            SAMRAI::hier::Box const& transformedSource,
                            SAMRAI::hier::Box const& destinationBox,
                            SAMRAI::hier::Box const& intersectionBox,
                            SAMRAI::hier::Transformation const& transformation,
                            ParticlesData const& source)
    {
        SAMRAI::hier::IntVector shift{transformation.getOffset()};
        SAMRAI::hier::Box intersectLocalSource{intersectionBox};

        transformation.inverseTransform(intersectLocalSource);

        intersectLocalSource.setLower(intersectLocalSource.lower() - sourceBox.lower());
        intersectLocalSource.setUpper(intersectLocalSource.upper() - sourceBox.lower());

        shift += sourceBox.lower();

        copy_(source, shift, intersectLocalSource);
    }


    void copy_(ParticlesData const& source, SAMRAI::hier::IntVector const& shiftToDestination,
               SAMRAI::hier::Box const& localSource)
    {
        for (auto const& particle : source.interior)
        {
            if (isInBox_(localSource, particle))
            {
                auto shiftedParticle{particle};

                shiftParticle_(shiftToDestination, shiftedParticle);

                if (isInBox_(interiorBox_, shiftedParticle))
                {
                    interior.push_back(std::move(shiftedParticle));
                }
                else
                {
                    ghost.push_back(std::move(shiftedParticle));
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


        localSource.setLower(localSource.lower() - sourceBox.lower());
        localSource.setUpper(localSource.upper() - sourceBox.lower());

        auto const& shift = sourceBox.lower();


        for (auto const& particle : interior)
        {
            if (isInBox_(localSource, particle))
            {
                auto shiftedParticle{particle};
                shiftParticle_(shift, shiftedParticle);
                buffer.push_back(shiftedParticle);
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
