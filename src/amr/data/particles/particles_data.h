#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_H
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_H

#include <numeric>
#include <stdexcept>

#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/MemoryUtilities.h>


#include "data/ions/ion_population/particle_pack.h"
#include "data/particles/particle.h"
#include "data/particles/particle_array.h"
#include "tools/amr_utils.h"

namespace PHARE
{
template<std::size_t dim>
static bool isInBox(SAMRAI::hier::Box const& box, Particle<dim> const& particle);

template<>
bool isInBox(SAMRAI::hier::Box const& box, Particle<1> const& particle)
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
inline bool isInBox(SAMRAI::hier::Box const& box, Particle<2> const& particle)
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
inline bool isInBox(SAMRAI::hier::Box const& box, Particle<3> const& particle)
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


/*template<>
bool isInBox(SAMRAI::hier::Box const& box, Particle<1> const& particle);

template<>
bool isInBox(SAMRAI::hier::Box const& box, Particle<1> const& particle);

template<>
bool isInBox(SAMRAI::hier::Box const& box, Particle<1> const& particle);*/


/** @brief ParticlesData is a concrete SAMRAI::hier::PatchData subclass to store Particle data
 *
 * This class encapsulates particle storage known by the module core, and by being derived
 * from PatchData is compatible with the SAMRAI data management system.
 *
 * A ParticlesData encapsulates **three** different particle arrays:
 *
 * - domainParticles : these particles are those for which iCell is within the physical domain of
 * the patch
 *
 * - ghostParticles: these particles are located within the ghost layer around the physical domain
 * of the patch. We call the "ghost layer" the layer of ghostCellWidth just outside the physical
 * domain of the patch, on borders that have neighbors patchs of the same level. Patch boundaries
 * that are at the level boundary are not ghost layers but coarseToFine boundaries. All the
 * particles in the ghost layer are exact clones of particles located on a neighbor patch of the
 * same level. The ghost particles are getting here when then exit the neighbor patch, and can enter
 * the patch.
 *
 * - coarseToFineParticles: these particles are located in a layer just passed the patch boundaries
 * that also are level boundaries. These particles are getting here when there is a particle
 * refinement from a coarser level
 *
 */
template<std::size_t dim>
class ParticlesData : public SAMRAI::hier::PatchData
{
public:
    ParticlesData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost)
        : SAMRAI::hier::PatchData::PatchData(box, ghost)
        , pack{&domainParticles, &ghostParticles, &coarseToFineParticles}
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
    virtual void copy(SAMRAI::hier::PatchData const& source) final
    {
        TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source);

        const ParticlesData* pSource = dynamic_cast<const ParticlesData*>(&source);
        if (pSource != nullptr)
        {
            SAMRAI::hier::Box const& sourceGhostBox = pSource->getGhostBox();
            SAMRAI::hier::Box const& myGhostBox     = getGhostBox();
            const SAMRAI::hier::Box intersectionBox{sourceGhostBox * myGhostBox};

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



    /**
     * @brief our copy2 will be called by a PatchData if a ParticleData was
     * given to be copied into another kind of PatchData. Here we chose that
     * copy2 throws unconditiionnally.
     */
    virtual void copy2(SAMRAI::hier::PatchData& destination) const final
    {
        throw std::runtime_error("Cannot cast");
    }



    /**
     * @brief copy with an overlap. Does the copy as the other overload but this time
     * the copy must account for the intersection with the boxes within the overlap
     * The copy is done between the source patch data and myself
     */
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
                    SAMRAI::hier::Box sourceGhostBox      = pSource->getGhostBox();
                    SAMRAI::hier::Box destinationGhostBox = this->getGhostBox();
                    SAMRAI::hier::Box intersectionBox{sourceGhostBox.getDim()};

                    if (isSameBlock(transformation))
                    {
                        if (offsetIsZero(transformation))
                        {
                            intersectionBox = overlapBox * sourceGhostBox * destinationGhostBox;

                            if (!intersectionBox.empty())
                            {
                                copy_(sourceGhostBox, destinationGhostBox, intersectionBox,
                                      *pSource);
                            }
                        }
                        else
                        {
                            SAMRAI::hier::Box shiftedSourceBox{sourceGhostBox};
                            transformation.transform(shiftedSourceBox);
                            intersectionBox = overlapBox * shiftedSourceBox * destinationGhostBox;

                            if (!intersectionBox.empty())
                            {
                                /* copy_(sourceGhostBox, destinationGhostBox, intersectionBox,
                                 *pSource);*/

                                copyWithTransform_(sourceGhostBox, intersectionBox, transformation,
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

                auto const& sourceBox = getGhostBox();

                SAMRAI::hier::Box transformedSource{sourceBox};
                transformation.transform(transformedSource);


                for (auto const& destinationBox : boxContainer)
                {
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

                auto myBox      = getBox();
                auto myGhostBox = getGhostBox();

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
                    // auto intersectLocalSource = AMRToLocal(intersect, getGhostBox());
                    // auto particleShift        = AMRToLocal(getGhostBox());

                    for (auto const& particle : particleArray)
                    {
                        // shift the particle to local index space
                        // and if it is in intersection, decide in which array to push it.
                        // auto shiftedParticle{particle};
                        // shiftParticle_(particleShift, shiftedParticle);

                        if (isInBox(intersect, particle))
                        {
                            if (isInBox(myBox, particle))
                            {
                                domainParticles.push_back(std::move(particle));
                            }
                            else
                            {
                                ghostParticles.push_back(std::move(particle));
                            }
                        }
                    } // end species loop
                }     // end box loop
            }         // end no rotation
        }             // end overlap not empty
    }



    ParticlesPack<ParticleArray<dim>>* getPointer() { return &pack; }



    // Core interface
    // these particles arrays are public because core module is free to use
    // them easily
    ParticleArray<dim> domainParticles;
    ParticleArray<dim> ghostParticles;

    ParticleArray<dim> coarseToFineParticles;

    ParticleArray<dim> coarseToFineParticlesOld;
    ParticleArray<dim> coarseToFineParticlesNew;

    ParticlesPack<ParticleArray<dim>> pack;



private:
    //! interiorLocalBox_ is the box, in local index space, that goes from the first to the last
    //! cell in our patch physical domain, i.e. "from dual physical start index to dual physical end
    //! index"
    SAMRAI::hier::Box interiorLocalBox_;




    void copy_(SAMRAI::hier::Box const& sourceGhostBox,
               SAMRAI::hier::Box const& destinationGhostBox,
               SAMRAI::hier::Box const& intersectionBox, ParticlesData const& sourceData)
    {
        std::array<decltype(sourceData.domainParticles) const*, 2> particlesArrays{
            &sourceData.domainParticles, &sourceData.ghostParticles};

        auto myDomainBox = this->getBox();


        // for each particles in the source ghost and domain particle arrays
        // we check if it is in the intersectionBox
        // if it is, is it in my domain box ?
        //      - if so, let's add it to my domain particle array
        //      - if not, let's add it to my ghost particle array
        for (auto const& sourceParticlesArray : particlesArrays)
        {
            for (auto const& particle : *sourceParticlesArray)
            {
                if (isInBox(intersectionBox, particle))
                {
                    if (isInBox(myDomainBox, particle))
                    {
                        domainParticles.push_back(particle);
                    }
                    else
                    {
                        ghostParticles.push_back(particle);
                    }
                }
            }
        }
    }




    void copyWithTransform_(SAMRAI::hier::Box const& sourceGhostBox,
                            SAMRAI::hier::Box const& intersectionBox,
                            SAMRAI::hier::Transformation const& transformation,
                            ParticlesData const& sourceData)
    {
        std::array<decltype(sourceData.domainParticles) const*, 2> particlesArrays{
            &sourceData.domainParticles, &sourceData.ghostParticles};

        auto myDomainBox = this->getBox();

        auto offset = transformation.getOffset();

        for (auto const& sourceParticlesArray : particlesArrays)
        {
            for (auto const& particle : *sourceParticlesArray)
            {
                // the particle is only copied if it is in the intersectionBox
                // but before its iCell must be shifted by the transformation offset

                auto newParticle{particle};
                for (auto iDir = 0; iDir < newParticle.iCell.size(); ++iDir)
                {
                    newParticle.iCell[iDir] += offset[iDir];
                }

                if (isInBox(intersectionBox, newParticle))
                {
                    // now we now the particle is in the intersection
                    // we need to know whether it is in the domain part of that
                    // intersection. If it is not, then it must be in the ghost part


                    if (isInBox(myDomainBox, newParticle))
                    {
                        domainParticles.push_back(newParticle);
                    }
                    else
                    {
                        ghostParticles.push_back(newParticle);
                    }
                }
            }
        }


        // SAMRAI::hier::Box localSourceSelectionBox = AMRToLocal(intersectionBox, sourceGhostBox);

        // we shift it back the box on top of source AMR indexes
        // transformation.inverseTransform(localSourceSelectionBox);

        // copy_(sourceData, particleShift, localSourceSelectionBox);
    }




    /**
     * @brief countNumberParticlesIn_ counts the number of interior particles that lie
     * within the boxes of an overlap
     */
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
            SAMRAI::hier::Box shiftedBox{box};
            SAMRAI::hier::Transformation const& transformation = overlap.getTransformation();
            transformation.transform(shiftedBox);
            SAMRAI::hier::Box intersectionBox{shiftedBox * getBox()};

            numberParticles += countNumberParticlesIn_(intersectionBox);
        }
        return numberParticles;
    }


    /**
     * @brief countNumberParticlesIn_ returns the number of interior particles within a given box
     *
     * the box given is in AMR index space so the function first needs to put it in
     * local indexing relative to the domain box
     */
    std::size_t countNumberParticlesIn_(SAMRAI::hier::Box const& box) const
    {
        std::size_t numberParticles{0};

        auto localSourceBox = AMRToLocal(box, getBox());

        for (auto const& particle : domainParticles)
        {
            if (isInBox(localSourceBox, particle))
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
        std::array<decltype(domainParticles) const*, 2> particlesArrays{&domainParticles,
                                                                        &ghostParticles};

        for (auto const& sourceParticlesArray : particlesArrays)
        {
            for (auto const& particle : *sourceParticlesArray)
            {
                auto shiftedParticle{particle};
                auto offset = transformation.getOffset();
                for (auto i = 0; i < dim; ++i)
                {
                    shiftedParticle.iCell[i] += offset[i];
                }
                if (isInBox(intersectionBox, shiftedParticle))
                {
                    buffer.push_back(shiftedParticle);
                }
            }
        }
    }
};


} // namespace PHARE

#endif
