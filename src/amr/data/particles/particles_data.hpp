#ifndef PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP
#define PHARE_SRC_AMR_DATA_PARTICLES_PARTICLES_DATA_HPP



#include "core/def.hpp" // IWYU pragma: keep
#include "core/logger.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_service.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"

#include "amr/samrai.hpp" // IWYU pragma: keep
#include "amr/utilities/box/amr_box.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include <SAMRAI/hier/BoxContainer.h>
#include <amr/data/particles/particles_variable_fill_pattern.hpp>


#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/RestartManager.h>
#include "SAMRAI/hier/Transformation.h"
#include <SAMRAI/tbox/MemoryUtilities.h>

#include <tuple>
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
    template<typename ParticleArray_>
    class ParticlesData : public SAMRAI::hier::PatchData
    {
        using This          = ParticlesData<ParticleArray_>;
        using Super         = SAMRAI::hier::PatchData;
        using ParticleArray = ParticleArray_;
        using SamBox        = SAMRAI::hier::Box;

        // using Particle_t          = typename ParticleArray::Particle_t;
        static constexpr auto dim = ParticleArray::dimension;
        // add one cell surrounding ghost box to map particles exiting the ghost layer
        static constexpr int ghostSafeMapLayer = 1;


        auto construct_particles(auto const& ghosts) const
        {
            if constexpr (any_in(ParticleArray::layout_mode, core::LayoutMode::AoSMapped))
                return ParticleArray{grow(phare_box_from<dim>(getGhostBox()), ghostSafeMapLayer)};
            else
                return make_particles<ParticleArray>(phare_box_from<dim>(getBox()), ghosts);
        }


        void validate_ghosts(auto const ghost)
        {
            core::for_N<dim>([&](auto i) {
                if (ghost[i] != ghost[0])
                    throw std::runtime_error("invalid");
            });
        }

    public:
        ParticlesData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost,
                      std::string const& name)
            : SAMRAI::hier::PatchData::PatchData(box, ghost)
            , domainParticles{construct_particles(ghost[0])}
            , patchGhostParticles{construct_particles(ghost[0])}
            , levelGhostParticles{construct_particles(ghost[0])}
            , levelGhostParticlesOld{construct_particles(ghost[0])}
            , levelGhostParticlesNew{construct_particles(ghost[0])}
            , pack{name,
                   &domainParticles,
                   &patchGhostParticles,
                   &levelGhostParticles,
                   &levelGhostParticlesOld,
                   &levelGhostParticlesNew}
            , interiorLocalBox_{AMRToLocal(box, this->getGhostBox())}
            , name_{name}
        {
            validate_ghosts(ghost);
        }


        auto& name() const { return name_; }

        ParticlesData()                                = delete;
        ParticlesData(ParticlesData const&)            = delete;
        ParticlesData(ParticlesData&&)                 = default;
        ParticlesData& operator=(ParticlesData const&) = delete;

        // SAMRAI interface
        void putToRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) const override
        {
            Super::putToRestart(restart_db);

            putParticlesToRestart(*restart_db, "domainParticles", domainParticles);
            putParticlesToRestart(*restart_db, "levelGhostParticles", levelGhostParticles);
            putParticlesToRestart(*restart_db, "levelGhostParticlesNew", levelGhostParticlesNew);
            putParticlesToRestart(*restart_db, "levelGhostParticlesOld", levelGhostParticlesOld);
        };


        void getFromRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) override
        {
            Super::getFromRestart(restart_db);

            getParticlesFromRestart(*restart_db, "domainParticles", domainParticles);
            getParticlesFromRestart(*restart_db, "levelGhostParticles", levelGhostParticles);
            getParticlesFromRestart(*restart_db, "levelGhostParticlesNew", levelGhostParticlesNew);
            getParticlesFromRestart(*restart_db, "levelGhostParticlesOld", levelGhostParticlesOld);
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

            SamBox const& sourceBox  = pSource.getBox();
            SamBox const& myGhostBox = getGhostBox();
            SamBox const intersectionBox{sourceBox * myGhostBox};

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
            return sizeof(std::size_t) + numberParticles * ParticleArray::size_of_particle();
        }


        std::size_t getCellOverlapDataStreamSize(SAMRAI::pdat::CellOverlap const& pOverlap) const
        {
            return sizeof(std::size_t)
                   + countNumberParticlesIn_(pOverlap) * ParticleArray::size_of_particle();
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
            using PackArray = core::ParticleArray<ParticleArray::options.with_layout(
                core::base_layout_type<ParticleArray>())>;

            if (pOverlap.isOverlapEmpty())
            {
                constexpr std::size_t zero = 0;
                stream << zero;
            }
            else
            {
                SAMRAI::hier::Transformation const& transformation = pOverlap.getTransformation();
                PackArray outBuffer; // make thread static
                pack_(pOverlap, transformation, outBuffer);
                stream << outBuffer.size();
                stream.growBufferAsNeeded();
                this->pack_(stream, outBuffer);
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


        template<typename ParticleArray_t>
        void pack_(SAMRAI::tbox::MessageStream& stream, ParticleArray_t const& outBuffer) const
        {
            if constexpr (any_in(ParticleArray_t::layout_mode, core::LayoutMode::SoA))
                std::apply(
                    [&](auto const&... container) {
                        ((stream.pack(container.data(), outBuffer.size())), ...);
                    },
                    outBuffer.as_tuple());
            else
                stream.pack(outBuffer.data(), outBuffer.size());
        }



        void unpack_from_ghost(SAMRAI::tbox::MessageStream& stream,
                               ParticlesDomainOverlap const& overlap);

        void unpack_cell_overlap(SAMRAI::tbox::MessageStream& stream,
                                 SAMRAI::pdat::CellOverlap const& pOverlap)
        {
            using UnpackArray = core::ParticleArray<ParticleArray::options.with_layout(
                core::base_layout_type<ParticleArray>())>;

            if (pOverlap.isOverlapEmpty())
                return;

            if (pOverlap.getTransformation().getRotation()
                != SAMRAI::hier::Transformation::NO_ROTATE)
                return;

            std::size_t numberParticles = 0;
            stream >> numberParticles;
            UnpackArray particleArray(numberParticles);
            unpack_(stream, particleArray);

            for (auto const& overlapBox : pOverlap.getDestinationBoxContainer())
            {
                auto const intersect = getGhostBox() * overlapBox;
                for (auto const& particle : particleArray)
                    if (isInBox(intersect, particle))
                        domainParticles.push_back(particle);
            }
            core::ParticleArrayService::template sync<0, core::ParticleType::Domain>(
                domainParticles);
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

        template<typename... Args>
        void unpack_from_ghost(Args&&... args);


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

        template<typename ParticleArray_t>
        void unpack_(SAMRAI::tbox::MessageStream& stream, ParticleArray_t& specie) const
        {
            using enum core::LayoutMode;
            if constexpr (any_in(ParticleArray_t::layout_mode, SoA /*, SoATS*/))
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
        SamBox interiorLocalBox_;


        std::string name_;

        void copy_(SamBox const& overlapBox, ParticlesData const& sourceData)
        {
            core::select_particles(sourceData.domainParticles, domainParticles,
                                   phare_box_from<dim>(overlapBox));
        }

        void copy_(SamBox const& overlapBox, ParticlesData const& sourceData,
                   SAMRAI::hier::Transformation const& transformation)
        {
            auto const offset = as_point<dim>(transformation);
            core::select_particles(sourceData.domainParticles, domainParticles,
                                   shift(phare_box_from<dim>(overlapBox), offset * -1), offset);
        }


        /**
         * @brief countNumberParticlesIn_ counts the number of particles that lie
         * within the boxes of an overlap. This function count both patchGhost and
         * domain particles since both could be streamed and we want an upperbound
         * on the number of bytes that could be streamed.
         */
        std::size_t countNumberParticlesIn_(SAMRAI::pdat::CellOverlap const& overlap) const
        {
            // throw std::runtime_error("This is never called!");
            // but if it is, below should work (maybe)

            PHARE_LOG_SCOPE(3, "ParticleData::countNumberParticlesIn_");

            if (overlap.isOverlapEmpty())
                return 0;

            SAMRAI::hier::Transformation const& transformation = overlap.getTransformation();
            std::size_t numberParticles                        = 0;
            for (auto const& overlapBox : overlap.getDestinationBoxContainer())
            {
                SamBox shiftedOverlapBox{overlapBox};
                transformation.inverseTransform(shiftedOverlapBox);
                numberParticles += core::count_particles(domainParticles,
                                                         phare_box_from<dim>(shiftedOverlapBox));
            }

            return numberParticles;
        }

        template<typename ParticleArray_t>
        void pack_(SAMRAI::pdat::CellOverlap const& overlap,
                   SAMRAI::hier::Transformation const& transformation,
                   ParticleArray_t& outBuffer) const
        {
            PHARE_LOG_SCOPE(3, "ParticleData::pack_");
            auto const offset = as_point<dim>(transformation);
            for (auto const& box : overlap.getDestinationBoxContainer())
                core::select_particles(domainParticles, outBuffer,
                                       shift(phare_box_from<dim>(box), offset * -1), offset);
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

    for (auto const& overlapBox : pOverlap.getDestinationBoxContainer())
        core::select_particles<core::ParticleType::PatchGhost>(
            src_particles, dst_particles, shift(phare_box_from<dim>(overlapBox), noffset), offset);
}




template<typename ParticleArray_t>
void ParticlesData<ParticleArray_t>::pack_from_ghost(SAMRAI::tbox::MessageStream& stream,
                                                     ParticlesDomainOverlap const& pOverlap) const
{
    using PackArray = core::ParticleArray<ParticleArray_t::options.with_layout(
        core::base_layout_type<ParticleArray_t>())>;

    PHARE_LOG_SCOPE(3, "ParticlesData::pack_from_ghost");

    if (pOverlap.isOverlapEmpty())
    {
        constexpr std::size_t zero = 0;
        stream << zero;
        return;
    }

    auto& src_particles = patchGhostParticles;
    auto const& offset  = as_point<dim>(pOverlap.getTransformation());
    auto const& noffset = offset * -1;

    PackArray outBuffer;
    for (auto const& overlapBox : pOverlap.getDestinationBoxContainer())
        core::select_particles<core::ParticleType::PatchGhost>(
            src_particles, outBuffer, shift(phare_box_from<dim>(overlapBox), noffset), offset);

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
    using PackArray = core::ParticleArray<ParticleArray_t::options.with_layout(
        core::base_layout_type<ParticleArray_t>())>;

    PHARE_LOG_SCOPE(3, "ParticlesData::unpack_from_ghost");

    std::size_t numberParticles = 0;
    stream >> numberParticles;
    PackArray particleArray(numberParticles);
    stream.unpack(particleArray.data(), numberParticles);

    for (auto const& p : particleArray)
        domainParticles.push_back(p);

    core::ParticleArrayService::template sync<0, core::ParticleType::Domain>(domainParticles);
}



} // namespace PHARE::amr



#endif
