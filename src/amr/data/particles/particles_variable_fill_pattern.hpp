#ifndef PHARE_SRC_AMR_PARTICLES_PARTICLES_VARIABLE_FILL_PATTERN_HPP
#define PHARE_SRC_AMR_PARTICLES_PARTICLES_VARIABLE_FILL_PATTERN_HPP

#include "core/logger.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "amr/utilities/box/amr_box.hpp"

#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/xfer/VariableFillPattern.h>

#include <vector>
#include <cassert>
#include <utility>
#include <stdexcept>
#include <algorithm>

namespace PHARE::amr
{


/** ParticlesDomainOverlap is used as a signal in particles_data.hpp
  that we are performing an export from the patch ghost layer
  of one patch, to the domain of adjacent patch which is the analogue
  of the original patch ghost layer
  */

class ParticlesDomainOverlap : public SAMRAI::pdat::CellOverlap
{
    using Super = SAMRAI::pdat::CellOverlap;

public:
    ParticlesDomainOverlap(SAMRAI::hier::BoxContainer const& boxes,
                           SAMRAI::hier::Transformation const& transformation)
        : Super{boxes, transformation}
    {
    }

    ~ParticlesDomainOverlap() = default;
};




/**
 * \brief VariableFillPattern that is used to grab particles leaving neighboring patches
 *
 * This pattern is only use to grab incoming particles and not in refinement operations.
 * Thus, only calculateOverlap is implemented. computeFillBoxesOverlap, only used in refinement
 * is not to be used.
 *
 * Leaving neighbor particles will be searched in a layer around neighbor patch domain.
 * Typically, because a particle should not travel more than one cell in one time step
 * this layer should be one cell wide.
 *
 * Here, we compute the overlap using the particle data geometry, which is the CellGeometry
 * This overlap will be the intersection of the source box with the destination ghost box.
 *
 * What we want is the opposite, we want to take leaving particles from the source GHOST box
 * that are found in the destination box.
 *
 * We thus grow the overlap by some amount to extend beyond the source box
 * and then we intersect with the destination box, to exclude the destination ghost layer.
 *
 * As explained above, the overlap could probably be grown by only one cell.
 * Here we use the particle ghost width defined in the GridLayout_t.
 */
template<typename GridLayout_t>
class ParticleDomainFromGhostFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    std::size_t constexpr static dim         = GridLayout_t::dimension;
    bool constexpr static overwrite_interior = false;

public:
    ParticleDomainFromGhostFillPattern() {}

    virtual ~ParticleDomainFromGhostFillPattern() {}

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& dst_geometry,
                     SAMRAI::hier::BoxGeometry const& src_geometry,
                     SAMRAI::hier::Box const& dst_patch_box, SAMRAI::hier::Box const& src_mask,
                     SAMRAI::hier::Box const& fill_box, bool const fn_overwrite_interior,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
        PHARE_LOG_SCOPE(3, "ParticleDomainFromGhostFillPattern::calculateOverlap");
#ifndef DEBUG_CHECK_DIM_ASSERTIONS
        NULL_USE(dst_patch_box);
#endif
        TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);

        auto basic_overlap = dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box,
                                                           overwrite_interior, transformation);

        auto& cell_overlap = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(*basic_overlap);

        SAMRAI::hier::BoxContainer boxes;
        for (auto const& box : cell_overlap.getDestinationBoxContainer())
        {
            auto const ghost_overlap
                = grow(phare_box_from<dim>(box), GridLayout_t::nbrParticleGhosts());
            auto const domain_overlap = ghost_overlap * phare_box_from<dim>(dst_patch_box);
            boxes.pushBack(samrai_box_from(*domain_overlap));
        }

        return std::make_shared<ParticlesDomainOverlap>(boxes, cell_overlap.getTransformation());
    }

    std::string const& getPatternName() const override { return s_name_id; }

private:
    ParticleDomainFromGhostFillPattern(ParticleDomainFromGhostFillPattern const&) = delete;
    ParticleDomainFromGhostFillPattern& operator=(ParticleDomainFromGhostFillPattern const&)
        = delete;

    static inline std::string const s_name_id = "BOX_GEOMETRY_FILL_PATTERN";

    SAMRAI::hier::IntVector const& getStencilWidth() override
    {
        TBOX_ERROR("getStencilWidth() should not be\n"
                   << "called.  This pattern creates overlaps based on\n"
                   << "the BoxGeometry objects and is not restricted to a\n"
                   << "specific stencil.\n");

        return SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension(1));
    }


    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(SAMRAI::hier::BoxContainer const& fill_boxes,
                            SAMRAI::hier::BoxContainer const& node_fill_boxes,
                            SAMRAI::hier::Box const& patch_box, SAMRAI::hier::Box const& data_box,
                            SAMRAI::hier::PatchDataFactory const& pdf) const override
    { throw std::runtime_error("no refinement supported or expected"); }
};

/**
 * \brief VariableFillPattern for level ghost particle filling that prevents duplicate fills.
 *
 * When two adjacent fine patches have overlapping ghost regions, the standard SAMRAI behavior
 * would fill both with coarse particles in their overlap — producing duplicates. This fill
 * pattern restricts each patch to only the portion of the level ghost region it exclusively
 * owns (no neighbor with lower flattened index claims it). The owned boxes must be set before
 * the RefineSchedule is built (i.e., before RefinerPool::registerLevel).
 */
template<typename GridLayoutT>
class ParticleLevelGhostVariableFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    static constexpr auto dim = GridLayoutT::dimension;

public:
    ParticleLevelGhostVariableFillPattern() {}

    virtual ~ParticleLevelGhostVariableFillPattern() {}

    void reset() { patchOwnedBoxes_.clear(); }

    void setOwnedBoxes(SAMRAI::hier::Box const& patchBox, SAMRAI::hier::BoxContainer const& owned)
    { patchOwnedBoxes_.emplace_back(patchBox, owned); }

    std::string const& getPatternName() const override { return s_name_id; }

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& dst_geometry,
                     SAMRAI::hier::BoxGeometry const& src_geometry,
                     SAMRAI::hier::Box const& /*dst_patch_box*/, SAMRAI::hier::Box const& src_mask,
                     SAMRAI::hier::Box const& fill_box, bool const overwrite_interior,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
        return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box, overwrite_interior,
                                             transformation);
    }

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(SAMRAI::hier::BoxContainer const& fill_boxes,
                            SAMRAI::hier::BoxContainer const& /*node_fill_boxes*/,
                            SAMRAI::hier::Box const& patch_box, SAMRAI::hier::Box const& data_box,
                            SAMRAI::hier::PatchDataFactory const& pdf) const override
    {
        SAMRAI::hier::Transformation transformation(
            SAMRAI::hier::IntVector::getZero(patch_box.getDim()));

        auto it = std::find_if(patchOwnedBoxes_.begin(), patchOwnedBoxes_.end(),
                               [&](auto const& p) { return p.first == patch_box; });

        SAMRAI::hier::BoxContainer overlap_boxes;
        if (it != patchOwnedBoxes_.end())
        {
            for (auto const& fillBox : fill_boxes)
                for (auto const& ownedBox : it->second)
                    if (auto const inter = fillBox * ownedBox; !inter.empty())
                        overlap_boxes.pushBack(inter);
        }
        else
        {
            overlap_boxes = fill_boxes;
            overlap_boxes.intersectBoxes(data_box);
        }

        return pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes, transformation);
    }

private:
    ParticleLevelGhostVariableFillPattern(ParticleLevelGhostVariableFillPattern const&) = delete;
    ParticleLevelGhostVariableFillPattern& operator=(ParticleLevelGhostVariableFillPattern const&)
        = delete;

    static inline std::string const s_name_id = "BOX_GEOMETRY_FILL_PATTERN";

    SAMRAI::hier::IntVector const& getStencilWidth() override
    {
        TBOX_ERROR("getStencilWidth() should not be called for "
                   "ParticleLevelGhostVariableFillPattern");
        return SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension(1));
    }



    std::vector<std::pair<SAMRAI::hier::Box, SAMRAI::hier::BoxContainer>> patchOwnedBoxes_{};
};


} // namespace PHARE::amr

#endif /* PHARE_SRC_AMR_PARTICLES_PARTICLES_VARIABLE_FILL_PATTERN_H */
