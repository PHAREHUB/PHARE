#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_H
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_H

#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/PatchData.h>

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
        throw std::runtime_error("Not Implemented Yet");
    }
    virtual void copy(SAMRAI::hier::PatchData const& source,
                      SAMRAI::hier::BoxOverlap const& overlap) final
    {
        throw std::runtime_error("Not Implemented Yet");
    }

    virtual void copy2(SAMRAI::hier::PatchData& destination,
                       SAMRAI::hier::BoxOverlap const& overlap) const final
    {
        throw std::runtime_error("Not Implemented Yet");
    }

    virtual bool canEstimateStreamSizeFromBox() const final
    {
        throw std::runtime_error("Not Implemented Yet");
    }

    virtual size_t getDataStreamSize(SAMRAI::hier::BoxOverlap const& overlap) const final
    {
        throw std::runtime_error("Not Implemented Yet");
    }
    virtual void packStream(SAMRAI::tbox::MessageStream& stream,
                            SAMRAI::hier::BoxOverlap const& overlap) const final
    {
        throw std::runtime_error("Not Implemented Yet");
    }

    virtual void unpackStream(SAMRAI::tbox::MessageStream& stream,
                              SAMRAI::hier::BoxOverlap const& overlap) final
    {
        throw std::runtime_error("Not Implemented Yet");
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
        SAMRAI::hier::IntVector shift{sourceBox.lower() - destinationBox.lower()};
        SAMRAI::hier::Box intersectLocalSource{intersectionBox};

        intersectLocalSource.setLower(intersectionBox.lower() - sourceBox.lower());
        intersectLocalSource.setUpper(intersectionBox.upper() - sourceBox.lower());

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


    bool isInBox_(SAMRAI::hier::Box const& box, Particle<dim> const& particle);


    void shiftParticle_(SAMRAI::hier::IntVector const& shift, Particle<dim>& particle);
};




template<>
bool ParticlesData<1>::isInBox_(SAMRAI::hier::Box const& box, Particle<1> const& particle)
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
bool ParticlesData<2>::isInBox_(SAMRAI::hier::Box const& box, Particle<2> const& particle)
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
bool ParticlesData<3>::isInBox_(SAMRAI::hier::Box const& box, Particle<3> const& particle)
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
void ParticlesData<1>::shiftParticle_(SAMRAI::hier::IntVector const& shift, Particle<1>& particle)
{
    particle.iCell[0] += shift[0];
}

template<>
void ParticlesData<2>::shiftParticle_(SAMRAI::hier::IntVector const& shift, Particle<2>& particle)
{
    particle.iCell[0] += shift[0];
    particle.iCell[1] += shift[1];
}
template<>
void ParticlesData<3>::shiftParticle_(SAMRAI::hier::IntVector const& shift, Particle<3>& particle)
{
    particle.iCell[0] += shift[0];
    particle.iCell[1] += shift[1];
    particle.iCell[2] += shift[2];
}

} // namespace PHARE

#endif
