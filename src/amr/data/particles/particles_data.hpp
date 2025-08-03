#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep



#include "core/def.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"

#include "amr/utilities/box/amr_box.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include <amr/data/particles/particles_variable_fill_pattern.hpp>


#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/RestartManager.h>
#include "SAMRAI/hier/Transformation.h"
#include <SAMRAI/tbox/MemoryUtilities.h>

#include <tuple>
#include <vector>
#include <cstddef>
#include <stdexcept>


namespace PHARE
{
namespace amr
{


    /** @brief ParticlesData is a concrete SAMRAI::hier::PatchData subclass
     * to store Particle data
     *
     * A ParticlesData encapsulates **three** different particle arrays:
     *
     * - domainParticles : these particles are those for which iCell
     *   is within the physical domain of the patch
     *
     * - patchGhostParticles: represents particles that left the patch domain and are
     *   physically located in the patch ghost layer of a patch.
     *
     * - levelGhostParticles: represent particles obtained from refinement and
     *   located in level ghost layer. These particles are to be pushed and injected
     *   in domain if they arrive in there.
     *
     *- levelGhostParticlesOld: same as levelGhostParticles but defined at previous
     *  next coarse time step. Used to deposit contribution of these particles
     *  to moments in level ghost nodes
     *
     *- levelGhostParticlesNew: same as levelGhostParticles but defined at next
     * coarser future time step. Used to deposit contribution of these particles
     * to moments in level ghost nodes
     *
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
        ParticlesData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost,
                      std::string const& name)
            : SAMRAI::hier::PatchData::PatchData(box, ghost)
            , domainParticles{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
            , patchGhostParticles{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
            , levelGhostParticles{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
            , levelGhostParticlesOld{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
            , levelGhostParticlesNew{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)}
            , pack{name,
                   &domainParticles,
                   &patchGhostParticles,
                   &levelGhostParticles,
                   &levelGhostParticlesOld,
                   &levelGhostParticlesNew}
            , interiorLocalBox_{AMRToLocal(box, this->getGhostBox())}
            , name_{name}
        {
        }


        auto& name() const { return name_; }

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
            putParticles("levelGhostParticles", levelGhostParticles);
            putParticles("levelGhostParticlesNew", levelGhostParticlesNew);
            putParticles("levelGhostParticlesOld", levelGhostParticlesOld);
        };


        void getFromRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) override
        {
            Super::getFromRestart(restart_db);

            using Packer = core::ParticlePacker<dim>;

            auto getParticles = [&](std::string const name, auto& particles) {
                std::array<bool, Packer::n_keys> keys_exist = core::generate(
                    [&](auto const& key) { return restart_db->keyExists(name + "_" + key); },
                    Packer::keys());

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



        template<typename... Args>
        void copy_from_ghost(Args&&... args);


        void copy_from_cell_overlap(ParticlesData const& pSource,
                                    SAMRAI::pdat::CellOverlap const& pOverlap)
        {
            SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
            SAMRAI::hier::BoxContainer const& boxList = pOverlap.getDestinationBoxContainer();
            for (auto const& overlapBox : boxList)
                copy_(overlapBox, pSource, transformation);
        }

        /**
         * @brief copy with an overlap given by SAMARAI.
         * At runtime we can deal with two kinds of overlaps:
         * - ParticlesDomainOverlap: means this copy is from a context when we're grabbing
         *   leaving domain particles from the neighbor patch, in the patchghost array.
         * - CellOverlap: means domain particles are copied as part of a refinement operation.
         */
        void copy(SAMRAI::hier::PatchData const& source,
                  SAMRAI::hier::BoxOverlap const& overlap) override
        {
            PHARE_LOG_SCOPE(3, "ParticlesData::copy with overlap");

            // casts throw on failure
            auto& pSource = dynamic_cast<ParticlesData const&>(source);

            if (auto particleOverlap = dynamic_cast<ParticlesDomainOverlap const*>(&overlap))
                copy_from_ghost(pSource, *particleOverlap);

            else if (auto pOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap))
                copy_from_cell_overlap(pSource, *pOverlap);

            else
                throw std::runtime_error("Unknown overlap type");
        }


        void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination,
                   [[maybe_unused]] SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            throw std::runtime_error("Cannot cast");
        }



        bool canEstimateStreamSizeFromBox() const override { return false; }

        std::size_t getOutGoingDataStreamSize(ParticlesDomainOverlap const& pOverlap) const
        {
            auto& transformation        = pOverlap.getTransformation();
            auto const& offset          = as_point<dim>(transformation);
            auto const& noffset         = offset * -1;
            std::size_t numberParticles = 0;
            for (auto const& overlapBox : pOverlap.getDestinationBoxContainer())
                numberParticles += patchGhostParticles.nbr_particles_in(
                    shift(phare_box_from<dim>(overlapBox), noffset));
            return sizeof(std::size_t) + numberParticles * sizeof(Particle_t);
        }


        std::size_t getCellOverlapDataStreamSize(SAMRAI::pdat::CellOverlap const& pOverlap) const
        {
            return sizeof(std::size_t) + countNumberParticlesIn_(pOverlap) * sizeof(Particle_t);
        }

        std::size_t getDataStreamSize(SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            if (auto particleOverlap = dynamic_cast<ParticlesDomainOverlap const*>(&overlap))
                return getOutGoingDataStreamSize(*particleOverlap);

            else if (auto pOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap))
                return getCellOverlapDataStreamSize(*pOverlap);

            else
                throw std::runtime_error("Unknown overlap type");
        }




        void pack_from_ghost(SAMRAI::tbox::MessageStream&, ParticlesDomainOverlap const&) const;

        void pack_from_cell_overlap(SAMRAI::tbox::MessageStream& stream,
                                    SAMRAI::pdat::CellOverlap const& pOverlap) const
        {
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
         * @brief packStream is the function that takes particles from our particles arrays
         * that lie in the boxes of the given overlap, and pack them to a stream.
         *
         * Streaming particles means that we have to take particles with iCell on a local source
         * index space , communicate them, and load them at destination with iCell on a
         * destination local index space. To do that we need to:
         *
         * 1- translate source iCell to source AMR index space
         * 2- Apply the offset to shift this AMR index on top of the destination cells
         * 3- pack and communicate particles
         * 4- move back iCell from the shifted AMR index space to the local destination index
         * space
         *
         * Note that step 2 could be done upon reception of the pack, we chose to do it before.
         *
         * As for copy(), we can have two kinds of overlaps:
         * - ParticlesDomainOverlap : for grabbing leaving domain particles
         * - CellOverlap : copy as part of refinement operations
         */
        void packStream(SAMRAI::tbox::MessageStream& stream,
                        SAMRAI::hier::BoxOverlap const& overlap) const override
        {
            PHARE_LOG_SCOPE(3, "ParticleData::packStream");

            if (auto particleOverlap = dynamic_cast<ParticlesDomainOverlap const*>(&overlap))
            {
                pack_from_ghost(stream, *particleOverlap);
            }
            else if (auto pOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap))
                pack_from_cell_overlap(stream, *pOverlap);
            else
                throw std::runtime_error("Unknown overlap type");
        }




        void unpack_from_ghost(SAMRAI::tbox::MessageStream& stream,
                               ParticlesDomainOverlap const& overlap);

        void unpack_cell_overlap(SAMRAI::tbox::MessageStream& stream,
                                 SAMRAI::pdat::CellOverlap const& pOverlap)
        {
            if (!pOverlap.isOverlapEmpty())
            {
                std::size_t numberParticles = 0;
                stream >> numberParticles;
                std::vector<Particle_t> particleArray(numberParticles);
                stream.unpack(particleArray.data(), numberParticles);


                SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
                if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
                {
                    SAMRAI::hier::BoxContainer const& overlapBoxes
                        = pOverlap.getDestinationBoxContainer();

                    for (auto const& overlapBox : overlapBoxes)
                    {
                        // note that we intersect the overlap box with the *ghost* box
                        // and not with the box although in the code, we never fill
                        // the patch ghost layer with particles (the level ghost layer
                        // is filled with particles but that is done in the refinement op).
                        // The reason for taking the ghost box YET putting the particles
                        // in the domain particle array is that SAMRAI may ask us to stream
                        // particles from a distant patch into a local temporary patch
                        // whost ghost box extends over the source data selection box.
                        // particles falling into our "ghost" layer here are thus not really
                        // ghost particles so they are just put in domain.
                        // Consistently, the ParticleRefineOperator will only look for
                        // particles to split from the domain particle array
                        //
                        // Note: see issue #1026 this intersection and check with isInBox
                        // may not be useful if particles all fall into the domain anyway
                        auto const intersect = getGhostBox() * overlapBox;

                        for (auto const& particle : particleArray)
                            if (isInBox(intersect, particle))
                                domainParticles.push_back(particle);

                    } // end box loop
                } // end no rotation
            } // end overlap not empty
        }

        /**
         * @brief unpackStream is the function that unpacks a stream of particles to our
         * domain particle array
         *
         * We get a stream and an overlap. The overlap contains boxes where to put particles and
         * transformation from source to destination AMR indexes.
         *
         * By convention chosen in patckStream, packed particles have their iCell in our AMR
         * index space since we are the destination.
         *
         * like for packStream, we can have two kinds of overlaps:
         * - ParticlesDomainOverlap : for unpacking leaving domain particles
         * - CellOverlap : unpacking as part of refinement operations
         *
         */
        void unpackStream(SAMRAI::tbox::MessageStream& stream,
                          SAMRAI::hier::BoxOverlap const& overlap) override
        {
            PHARE_LOG_SCOPE(3, "ParticleData::unpackStream");

            if (auto* particleOverlap = dynamic_cast<ParticlesDomainOverlap const*>(&overlap))
                unpack_from_ghost(stream, *particleOverlap);

            else if (auto const* pOverlap
                     = dynamic_cast<SAMRAI::pdat::CellOverlap const*>(&overlap))
                unpack_cell_overlap(stream, *pOverlap);

            else
                throw std::runtime_error("Unknown overlap type");
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
        //! interiorLocalBox_ is the box, in local index space, that goes from the first to the
        //! last cell in our patch physical domain, i.e. "from dual physical start index to dual
        //! physical end index"
        SAMRAI::hier::Box interiorLocalBox_;
        std::string name_;

        void copy_(SAMRAI::hier::Box const& overlapBox, ParticlesData const& sourceData)
        {
            auto myDomainBox         = this->getBox();
            auto& srcDomainParticles = sourceData.domainParticles;

            PHARE_LOG_START(3, "ParticleData::copy_ DomainToDomain");

            // first copy particles that fall into our domain array
            // they can come from the source domain or patch ghost
            auto const destBox = myDomainBox * overlapBox;
            auto new_size      = domainParticles.size();

            if (!destBox.empty())
            {
                auto const destBox_p = phare_box_from<dim>(destBox);
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

            new_size = domainParticles.size();
            for (auto& selectionBox : ghostLayerBoxes)
            {
                if (!selectionBox.empty())
                {
                    auto selectionBox_p = phare_box_from<dim>(selectionBox);
                    new_size += srcDomainParticles.nbr_particles_in(selectionBox_p);
                }
            }
            if (domainParticles.capacity() < new_size)
                domainParticles.reserve(new_size);


            for (auto const& selectionBox : ghostLayerBoxes)
            {
                if (!selectionBox.empty())
                {
                    auto selectionBox_p = phare_box_from<dim>(selectionBox);
                    srcDomainParticles.export_particles(selectionBox_p, domainParticles);
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

            new_size = domainParticles.size();
            for (auto& selectionBox : ghostLayerBoxes)
            {
                if (!selectionBox.empty())
                {
                    transformation.inverseTransform(selectionBox);
                    auto selectionBox_p = phare_box_from<dim>(selectionBox);
                    new_size += srcDomainParticles.nbr_particles_in(selectionBox_p);
                }
            }
            if (domainParticles.capacity() < new_size)
                domainParticles.reserve(new_size);


            // ghostLayer boxes already have been inverse transformed
            // in previous loop, not to do again...
            for (auto const& selectionBox : ghostLayerBoxes)
            {
                if (!selectionBox.empty())
                {
                    auto selectionBox_p = phare_box_from<dim>(selectionBox);
                    srcDomainParticles.export_particles(selectionBox_p, domainParticles, offseter);
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


namespace PHARE::amr
{

template<typename ParticleArray_t>
template<typename... Args>
void ParticlesData<ParticleArray_t>::copy_from_ghost(Args&&... args)
{
    PHARE_LOG_SCOPE(3, "ParticlesData::copy_from_ghost");

    auto&& [pSource, pOverlap] = std::forward_as_tuple(args...);
    auto& src_particles        = pSource.patchGhostParticles;
    auto& dst_particles        = domainParticles;
    auto const& offset         = as_point<dim>(pOverlap.getTransformation());
    auto const& noffset        = offset * -1;

    auto const offsetToDest = [&](auto const& particle) {
        auto shiftedParticle{particle};
        for (std::size_t idir = 0; idir < dim; ++idir)
            shiftedParticle.iCell[idir] += offset[idir];
        return shiftedParticle;
    };
    // we shift the overlap box to the our array index space since it is given
    // in the destinaton index space.
    for (auto const& overlapBox : pOverlap.getDestinationBoxContainer())
        src_particles.export_particles(shift(phare_box_from<dim>(overlapBox), noffset),
                                       dst_particles, offsetToDest);
}



template<typename ParticleArray_t>
void ParticlesData<ParticleArray_t>::pack_from_ghost(SAMRAI::tbox::MessageStream& stream,
                                                     ParticlesDomainOverlap const& pOverlap) const
{
    PHARE_LOG_SCOPE(3, "ParticlesData::pack_from_ghost");

    if (pOverlap.isOverlapEmpty())
    {
        constexpr std::size_t zero = 0;
        stream << zero;
        return;
    }

    std::vector<Particle_t> outBuffer;
    auto& src_particles = patchGhostParticles;
    auto const& offset  = as_point<dim>(pOverlap.getTransformation());
    auto const& noffset = offset * -1;

    auto const offsetToDest = [&](auto const& particle) {
        auto shiftedParticle{particle};
        for (std::size_t idir = 0; idir < dim; ++idir)
            shiftedParticle.iCell[idir] += offset[idir];
        return shiftedParticle;
    };

    // we shift the overlap box to the our array index space since it is given
    // in the destinaton index space.
    for (auto const& overlapBox : pOverlap.getDestinationBoxContainer())
        src_particles.export_particles(shift(phare_box_from<dim>(overlapBox), noffset), outBuffer,
                                       offsetToDest);

    stream << outBuffer.size();
    stream.growBufferAsNeeded();
    stream.pack(outBuffer.data(), outBuffer.size());
}

// The overlap is not needed here as the pack selects only from the desired overlap
//  and the transform if applicable is performed during packing
template<typename ParticleArray_t>
void ParticlesData<ParticleArray_t>::unpack_from_ghost(SAMRAI::tbox::MessageStream& stream,
                                                       ParticlesDomainOverlap const& /*pOverlap*/)
{
    PHARE_LOG_SCOPE(3, "ParticlesData::unpack_from_ghost");

    std::size_t numberParticles = 0;
    stream >> numberParticles;
    std::vector<Particle_t> particleArray(numberParticles);
    stream.unpack(particleArray.data(), numberParticles);

    domainParticles.reserve(domainParticles.size() + numberParticles);
    // we disregard the overlap boxes in this function
    // contrary to unpack_cell_overlap.
    // the reason is that we only get here when we're unpacking
    // particles that are leaving neighbor domain into and so they
    // must be in the domain box, no need to check.
    for (auto const& p : particleArray)
        domainParticles.push_back(p);
}

} // namespace PHARE::amr

#endif
