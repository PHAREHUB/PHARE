#ifndef PHARE_AMR_MAGNETIC_PATCH_STRATEGY_BASE_HPP
#define PHARE_AMR_MAGNETIC_PATCH_STRATEGY_BASE_HPP

#include "coarse_cell_round_out.hpp"

#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <cstddef>
#include <stdexcept>
#include <string>

namespace PHARE::amr
{
/**
 * @brief Non-templated polymorphic base for the magnetic RefinePatchStrategy family.
 *
 * Declares the small contract every magnetic postprocess strategy needs (registerIDs) and
 * provides the boilerplate they all share:
 *   - the no-op SAMRAI overrides (setPhysicalBoundaryConditions, preprocessRefine) — the magnetic
 *     path has no physical-boundary or pre-refine work of its own;
 *   - getRefineOpStencilWidth, always 0: every magnetic postprocess (the Tóth-Roe base term, the
 *     ADPT touch-up) reconstructs an interior face from *fine* faces of the coarse cell it
 *     belongs to — faces the stage-1 gather has already written — and reads no coarse data at
 *     all. SAMRAI provisions max(this, every registered refine operator), and the operators
 *     report their own coarse reach (1 at order 2), so 0 here narrows nothing;
 *   - reconstructionRegion, the whole-coarse-cell round-out every postprocessRefine must apply to
 *     its fill box before reconstructing (see coarse_cell_round_out.hpp for the invariant this
 *     enforces and why SAMRAI's fill boxes do not satisfy it on their own).
 */
class MagneticPatchStrategyBase : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    virtual void registerIDs(int b_id) = 0;

    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch&, double const,
                                       SAMRAI::hier::IntVector const&) override
    {
    }

    void preprocessRefine(SAMRAI::hier::Patch&, SAMRAI::hier::Patch const&,
                          SAMRAI::hier::Box const&, SAMRAI::hier::IntVector const&) override
    {
    }

    SAMRAI::hier::IntVector
    getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector::getZero(dim);
    }

    /**
     * @brief The region a magnetic interior-face reconstruction over the fill box `fine_box` must
     * run over.
     *
     * Every magnetic postprocess (Tóth-Roe base term, ADPT touch-up) reaches exactly one coarse
     * cell: it reads only fine faces bounding the coarse cell it is reconstructing, and every one
     * of those is a shared face the stage-1 gather fills. So the reconstruction is self-consistent
     * iff the region it runs over is a union of whole coarse cells (see coarse_cell_round_out.hpp).
     *
     * SAMRAI's fill boxes are not. A recursive schedule fills a coarse-interpolation temporary plus
     * a ring of d_max_stencil_width = 1 cell, and one cell is always *half* a coarse cell, whatever
     * the temporary's own parity — so the shared faces the reconstruction needs on the far side of
     * a cut cell were never written.
     *
     * Rounding out grows each edge by at most one index, so the region stays inside the allocation
     * iff field_ghost_width >= fill_ring + 1. That holds in every configuration we support — the
     * ring is 1 (the composite refiner is order 2 only), hybrid has 2 ghosts and MHD 6 — so the
     * clip below never actually bites. If it ever did it would cut a coarse cell back in half, and
     * a half-covered cell has no well-posed reconstruction: some of its inputs are faces nothing
     * ever wrote. That is a configuration we do not handle, so say so rather than silently leaving
     * NaN faces behind for something downstream to sample.
     */
    template<std::size_t dimension>
    static SAMRAI::hier::Box reconstructionRegion(SAMRAI::hier::Box const& fine_box,
                                                  SAMRAI::hier::Box const& ghost_box)
    {
        auto const region = roundCellBoxOutToCoarseCells<dimension>(fine_box) * ghost_box;

        if (!isWholeCoarseCells<dimension>(region))
            throw std::runtime_error(
                "magnetic prolongation region is not a union of whole coarse cells: fill box "
                + to_str(fine_box) + ", rounded out and clipped to the allocated "
                + to_str(ghost_box) + ", gives " + to_str(region)
                + ", which half-covers a coarse cell. The field ghost width must be at least the "
                  "fill ring plus one.");

        return region;
    }

    ~MagneticPatchStrategyBase() override = default;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_MAGNETIC_PATCH_STRATEGY_BASE_HPP
