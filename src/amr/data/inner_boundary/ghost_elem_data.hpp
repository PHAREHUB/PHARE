#ifndef PHARE_AMR_DATA_INNER_BOUNDARY_GHOST_ELEM_DATA_HPP
#define PHARE_AMR_DATA_INNER_BOUNDARY_GHOST_ELEM_DATA_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/inner_boundary/ghost_elem_pack.hpp"
#include "core/inner_boundary/inner_boundary_defs.hpp"

#include "amr/samrai.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/tbox/MessageStream.h>

#include <cstddef>
#include <string>

namespace PHARE::amr
{

/**
 * @brief PatchData holding the per-patch precomputed ghost-element vectors
 * for one inner boundary.
 *
 * The actual storage is a `std::array<std::vector<GhostElemData<dim>>, 2^dim>`.
 * Each vector is purely patch-local (filled by the classifier at level init,
 * consumed by the BC applier each time step) — no halo exchange, no refinement
 * communication, no restart persistence — so all SAMRAI communication virtuals
 * are no-ops.
 */
template<std::size_t dim>
class GhostElemPatchData : public SAMRAI::hier::PatchData
{
    using Super = SAMRAI::hier::PatchData;

public:
    using pack_type = core::GhostElemPack<dim>;
    using array_type = typename pack_type::ghost_elem_array_type;

    GhostElemPatchData(SAMRAI::hier::Box const& box, SAMRAI::hier::IntVector const& ghost,
                       std::string const& name)
        : Super::PatchData(box, ghost)
        , data_{}
        , pack_{name}
    {
        pack_._data = &data_;
    }

    GhostElemPatchData()                                     = delete;
    GhostElemPatchData(GhostElemPatchData const&)            = delete;
    GhostElemPatchData(GhostElemPatchData&&)                 = default;
    GhostElemPatchData& operator=(GhostElemPatchData const&) = delete;

    auto& name() const { return pack_._name; }

    pack_type* getPointer() { return &pack_; }

    // ---- SAMRAI PatchData interface ---------------------------------------

    void copy(SAMRAI::hier::PatchData const& /*source*/) override {}
    void copy2(SAMRAI::hier::PatchData& /*destination*/) const override {}
    void copy(SAMRAI::hier::PatchData const& /*source*/,
              SAMRAI::hier::BoxOverlap const& /*overlap*/) override
    {
    }
    void copy2(SAMRAI::hier::PatchData& /*destination*/,
               SAMRAI::hier::BoxOverlap const& /*overlap*/) const override
    {
    }

    bool canEstimateStreamSizeFromBox() const override { return true; }

    std::size_t getDataStreamSize(SAMRAI::hier::BoxOverlap const& /*overlap*/) const override
    {
        return 0;
    }

    void packStream(SAMRAI::tbox::MessageStream& /*stream*/,
                    SAMRAI::hier::BoxOverlap const& /*overlap*/) const override
    {
    }

    void unpackStream(SAMRAI::tbox::MessageStream& /*stream*/,
                      SAMRAI::hier::BoxOverlap const& /*overlap*/) override
    {
    }

    void putToRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) const override
    {
        Super::putToRestart(restart_db);
    }

    void getFromRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) override
    {
        Super::getFromRestart(restart_db);
    }

private:
    array_type data_;
    pack_type  pack_;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_DATA_INNER_BOUNDARY_GHOST_ELEM_DATA_HPP
