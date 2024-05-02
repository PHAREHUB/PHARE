#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP

#include <iterator>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/MemoryUtilities.h>
#include <SAMRAI/tbox/RestartManager.h>
#include "SAMRAI/hier/Transformation.h"


#include "core/def.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/utilities/box/amr_box.hpp"
#include "core/utilities/point/point.hpp"

#include "core/logger.hpp"



namespace PHARE
{
namespace amr
{


    template<typename Particle>
    NO_DISCARD inline bool isInBox(SAMRAI::hier::Box const& box, Particle const& particle)
    {
        constexpr auto dim = Particle::dimension;

        auto const& iCell = particle.iCell;

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
        using Super = SAMRAI::hier::PatchData;

        using Particle_t          = typename ParticleArray::Particle_t;
        static constexpr auto dim = ParticleArray::dimension;
        // add one cell surrounding ghost box to map particles exiting the ghost layer
        static constexpr int ghostSafeMapLayer = 1;

    public:
        ParticlesData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost)
            : SAMRAI::hier::PatchData::PatchData(box, ghost)
            , domainParticles{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
            , patchGhostParticles{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
            , levelGhostParticles{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
            , levelGhostParticlesOld{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
            , levelGhostParticlesNew{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
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

        void putToRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) const override
        {
            Super::putToRestart(restart_db);

            using Packer = core::ParticlePacker<dim>;

            auto putParticles = [&](std::string name, auto& particles) {
                // SAMRAI errors on writing 0 size arrays
                if (particles.size() == 0)
                    return;

                particles.sortMapping();

                Packer packer(particles);
                core::ContiguousParticles<dim> soa{particles.size()};
                packer.pack(soa);

                std::size_t part_idx = 0;
                core::apply(soa.as_tuple(), [&](auto const& arg) {
                    restart_db->putVector(name + "_" + packer.keys()[part_idx++], arg);
                });
            };

            putParticles("domainParticles", domainParticles);
            putParticles("patchGhostParticles", patchGhostParticles);
            putParticles("levelGhostParticles", levelGhostParticles);
            putParticles("levelGhostParticlesNew", levelGhostParticlesNew);
            putParticles("levelGhostParticlesOld", levelGhostParticlesOld);
        };


        void getFromRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) override
        {
            Super::getFromRestart(restart_db);

            using Packer = core::ParticlePacker<dim>;

            auto getParticles = [&](std::string const name, auto& particles) {
                auto const keys_exist = core::generate(
                    [&](auto const& key) { return restart_db->keyExists(name + "_" + key); },
                    Packer::keys());

                assert(keys_exist.size() > 0);

                bool all  = core::all(keys_exist);
                bool none = core::none(keys_exist);
                if (!(all or none))
                    throw std::runtime_error("ParticlesData::getFromRestart has been given an "
                                             "invalid input file, inconsistent state detected");

                if (none) // can't read what doesn't exist
                    return;

                auto n_particles
                    = restart_db->getArraySize(name + "_" + Packer::arbitrarySingleValueKey());
                core::ContiguousParticles<dim> soa{n_particles};

                {
                    std::size_t part_idx = 0;
                    core::apply(soa.as_tuple(), [&](auto& arg) {
                        restart_db->getVector(name + "_" + Packer::keys()[part_idx++], arg);
                    });
                }

                assert(particles.size() == 0);
                particles.reserve(n_particles);
                for (std::size_t i = 0; i < n_particles; ++i)
                    particles.push_back(soa.copy(i));
            };

            getParticles("domainParticles", domainParticles);
            getParticles("patchGhostParticles", patchGhostParticles);
            getParticles("levelGhostParticles", levelGhostParticles);
            getParticles("levelGhostParticlesNew", levelGhostParticlesNew);
            getParticles("levelGhostParticlesOld", levelGhostParticlesOld);
        }


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
            PHARE_LOG_SCOPE(3, "ParticlesData::copy");

            TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source);

            // throws if fails
            auto& pSource = dynamic_cast<ParticlesData const&>(source);

            SAMRAI::hier::Box const& sourceBox  = pSource.getBox();
            SAMRAI::hier::Box const& myGhostBox = getGhostBox();
            const SAMRAI::hier::Box intersectionBox{sourceBox * myGhostBox};

            if (!intersectionBox.empty())
            {
                copy_(intersectionBox, pSource);
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
            PHARE_LOG_SCOPE(3, "ParticlesData::copy with overlap");

            // casts throw on failure
            auto& pSource  = dynamic_cast<ParticlesData const&>(source);
            auto& pOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap);

            SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
            SAMRAI::hier::BoxContainer const& boxList = pOverlap.getDestinationBoxContainer();
            for (auto const& overlapBox : boxList)
            {
                copy_(overlapBox, pSource, transformation);
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

            return countNumberParticlesIn_(pOverlap) * sizeof(Particle_t);
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
            PHARE_LOG_SCOPE(3, "ParticleData::packStream");

            auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};

            std::vector<Particle_t> outBuffer;

            if (pOverlap.isOverlapEmpty())
            {
                constexpr std::size_t zero = 0;
                stream << zero;
            }
            else
            {
                SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
                SAMRAI::hier::BoxContainer const& boxContainer
                    = pOverlap.getDestinationBoxContainer();
                pack_(pOverlap, transformation, outBuffer);
                stream << outBuffer.size();
                stream.growBufferAsNeeded();
                stream.pack(outBuffer.data(), outBuffer.size());
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
        void unpackStream(SAMRAI::tbox::MessageStream& stream,
                          SAMRAI::hier::BoxOverlap const& overlap) override
        {
            PHARE_LOG_SCOPE(3, "ParticleData::unpackStream");

            auto const& pOverlap{dynamic_cast<SAMRAI::pdat::CellOverlap const&>(overlap)};

            if (!pOverlap.isOverlapEmpty())
            {
                // unpack particles into a particle array
                std::size_t numberParticles = 0;
                stream >> numberParticles;
                std::vector<Particle_t> particleArray(numberParticles);
                stream.unpack(particleArray.data(), numberParticles);

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

                        for (auto const& particle : particleArray)
                        {
                            if (isInBox(intersect, particle))
                            {
                                if (isInBox(myBox, particle))
                                {
                                    domainParticles.push_back(particle);
                                }
                                else
                                {
                                    patchGhostParticles.push_back(particle);
                                }
                            }
                        } // end species loop
                    }     // end box loop
                }         // end no rotation
            }             // end overlap not empty
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

        void copy_(SAMRAI::hier::Box const& overlapBox, ParticlesData const& sourceData)
        {
            auto myDomainBox         = this->getBox();
            auto& srcDomainParticles = sourceData.domainParticles;

            PHARE_LOG_START(3, "ParticleData::copy_ DomainToDomain");

            // first copy particles that fall into our domain array
            // they can come from the source domain or patch ghost
            auto destBox  = myDomainBox * overlapBox;
            auto new_size = domainParticles.size();

            if (!destBox.empty())
            {
                auto destBox_p = phare_box_from<dim>(destBox);
                new_size += srcDomainParticles.nbr_particles_in(destBox_p);
                if (domainParticles.capacity() < new_size)
                    domainParticles.reserve(new_size);

                srcDomainParticles.export_particles(destBox_p, domainParticles);
            }

            PHARE_LOG_START(3, "ParticlesData::copy_ DomainToGhosts");
            // Now copy particles from the source domain that fall into
            // our ghost layer. The ghost layer is the result of removing the domain box
            // from the intersection box.
            SAMRAI::hier::BoxContainer ghostLayerBoxes{};
            ghostLayerBoxes.removeIntersections(overlapBox, myDomainBox);

            new_size = patchGhostParticles.size();
            for (auto& selectionBox : ghostLayerBoxes)
            {
                if (!selectionBox.empty())
                {
                    auto selectionBox_p = phare_box_from<dim>(selectionBox);
                    new_size += srcDomainParticles.nbr_particles_in(selectionBox_p);
                }
            }
            if (patchGhostParticles.capacity() < new_size)
                patchGhostParticles.reserve(new_size);


            for (auto const& selectionBox : ghostLayerBoxes)
            {
                if (!selectionBox.empty())
                {
                    auto selectionBox_p = phare_box_from<dim>(selectionBox);
                    srcDomainParticles.export_particles(selectionBox_p, patchGhostParticles);
                }
            }
            PHARE_LOG_STOP(3, "ParticlesData::copy_ DomainToGhosts");
        }

        void copy_(SAMRAI::hier::Box const& overlapBox, ParticlesData const& sourceData,
                   SAMRAI::hier::Transformation const& transformation)
        {
            auto myDomainBox         = this->getBox();
            auto& srcDomainParticles = sourceData.domainParticles;

            PHARE_LOG_START(3, "ParticleData::copy_ (transform)");

            // first copy particles that fall into our domain array
            // they can come from the source domain or patch ghost
            auto destBox  = myDomainBox * overlapBox;
            auto new_size = domainParticles.size();
            auto offset   = transformation.getOffset();
            auto offseter = [&](auto const& particle) {
                // we make a copy because we do not want to
                // shift the original particle...
                auto shiftedParticle{particle};
                for (std::size_t idir = 0; idir < dim; ++idir)
                {
                    shiftedParticle.iCell[idir] += offset[idir];
                }
                return shiftedParticle;
            };

            PHARE_LOG_START(3, "DomainToDomain (transform)");
            if (!destBox.empty())
            {
                // we cannot select particles from the intersectDomain box
                // right away. The reason is that the transformation may have
                // a non-zero offset and particle iCells from the source are in
                // the source index space, not in the destination index space
                // therefore we need to first modify the destination box to
                // be in the source index space
                // this is done by applying the INVERSE transformation
                // since a *transformation* is from source to destination.

                transformation.inverseTransform(destBox);
                auto destBox_p = phare_box_from<dim>(destBox);
                new_size += srcDomainParticles.nbr_particles_in(destBox_p);

                if (domainParticles.capacity() < new_size)
                    domainParticles.reserve(new_size);
                srcDomainParticles.export_particles(destBox_p, domainParticles, offseter);
            }
            PHARE_LOG_STOP(3, "DomainToDomain (transform)");



            PHARE_LOG_START(3, "DomainToGhosts (transform)");
            // Now copy particles from the source domain and patchghost that fall into
            // our ghost layer. The ghost layer is the result of removing the domain box
            // from the intersection box.
            SAMRAI::hier::BoxContainer ghostLayerBoxes{};
            ghostLayerBoxes.removeIntersections(overlapBox, myDomainBox);

            new_size = patchGhostParticles.size();
            for (auto& selectionBox : ghostLayerBoxes)
            {
                if (!selectionBox.empty())
                {
                    transformation.inverseTransform(selectionBox);
                    auto selectionBox_p = phare_box_from<dim>(selectionBox);
                    new_size += srcDomainParticles.nbr_particles_in(selectionBox_p);
                }
            }
            if (patchGhostParticles.capacity() < new_size)
                patchGhostParticles.reserve(new_size);


            // ghostLayer boxes already have been inverse transformed
            // in previous loop, not to do again...
            for (auto const& selectionBox : ghostLayerBoxes)
            {
                if (!selectionBox.empty())
                {
                    auto selectionBox_p = phare_box_from<dim>(selectionBox);
                    srcDomainParticles.export_particles(selectionBox_p, patchGhostParticles,
                                                        offseter);
                }
            }

            PHARE_LOG_STOP(3, "DomainToGhosts (transform)");
            PHARE_LOG_STOP(3, "ParticleData::copy_ (transform)");
        }




        /**
         * @brief countNumberParticlesIn_ counts the number of particles that lie
         * within the boxes of an overlap. This function count both patchGhost and
         * domain particles since both could be streamed and we want an upperbound
         * on the number of bytes that could be streamed.
         */
        std::size_t countNumberParticlesIn_(SAMRAI::pdat::CellOverlap const& overlap) const
        {
            PHARE_LOG_SCOPE(3, "ParticleData::countNumberParticlesIn_");
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
                SAMRAI::hier::Box shiftedOverlapBox{overlapBox};
                SAMRAI::hier::Transformation const& transformation = overlap.getTransformation();
                transformation.inverseTransform(shiftedOverlapBox);
                auto shiftedOverlapBox_p = phare_box_from<dim>(shiftedOverlapBox);
                numberParticles += domainParticles.nbr_particles_in(shiftedOverlapBox_p);
            }
            return numberParticles;
        }



        void pack_(SAMRAI::pdat::CellOverlap const& overlap,
                   SAMRAI::hier::Transformation const& transformation,
                   std::vector<Particle_t>& outBuffer) const
        {
            PHARE_LOG_SCOPE(3, "ParticleData::pack_");
            // we want to put particles from our domain and patchghost arrays
            // that fall into the intersection box Note that the overlap boxes
            // are not in the same index space as our particles.  the
            // transformation offset goes from OUR index space to the
            // destination space.  Therefore we need to inverse transform the
            // overlap box into our index space, intersect each of them with
            // our ghost box and put export them with the transformation offset
            auto overlapBoxes = overlap.getDestinationBoxContainer();
            auto offset       = transformation.getOffset();
            std::size_t size  = 0;
            auto offseter     = [&](auto const& particle) {
                auto shiftedParticle{particle};
                for (std::size_t idir = 0; idir < dim; ++idir)
                {
                    shiftedParticle.iCell[idir] += offset[idir];
                }
                return shiftedParticle;
            };
            for (auto const& box : overlapBoxes)
            {
                auto toTakeFrom{box};
                transformation.inverseTransform(toTakeFrom);
                auto toTakeFrom_p = phare_box_from<dim>(toTakeFrom);
                size += domainParticles.nbr_particles_in(toTakeFrom_p);
            }
            outBuffer.reserve(size);
            for (auto const& box : overlapBoxes)
            {
                auto toTakeFrom{box};
                transformation.inverseTransform(toTakeFrom);
                auto toTakeFrom_p = phare_box_from<dim>(toTakeFrom);
                domainParticles.export_particles(toTakeFrom_p, outBuffer, offseter);
            }
        }
    };
} // namespace amr

} // namespace PHARE

#endif
