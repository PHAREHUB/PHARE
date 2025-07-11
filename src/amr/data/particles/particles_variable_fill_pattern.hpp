#ifndef PHARE_SRC_AMR_PARTICLES_PARTICLES_VARIABLE_FILL_PATTERN_HPP
#define PHARE_SRC_AMR_PARTICLES_PARTICLES_VARIABLE_FILL_PATTERN_HPP

#include "core/def/phare_mpi.hpp"
#include <amr/utilities/box/amr_box.hpp>

#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/hier/RefineOperator.h>
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <cassert>
#include <stdexcept>

namespace PHARE::amr
{

// We use this class as a signal in particles_data.hpp
//  that we are performing an export from the patch ghost layer
//  of one patch, to the domain of adjacent patch which is the analogue
//  of the original patch ghost layer
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
    {
        PHARE_LOG_SCOPE(2, "ParticleDomainFromGhostFillPattern::computeFillBoxesOverlap");

        throw std::runtime_error("no refinement supported or expected");
    }
};

} // namespace PHARE::amr

#endif /* PHARE_SRC_AMR_PARTICLES_PARTICLES_VARIABLE_FILL_PATTERN_H */
