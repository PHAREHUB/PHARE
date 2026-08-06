#ifndef PHARE_MAGNETIC_COMPOSITE_REFINER_HPP
#define PHARE_MAGNETIC_COMPOSITE_REFINER_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "field_refiner_kernel.hpp"
#include "composite_field_refiner.hpp"

#include <memory>
#include <stdexcept>


namespace PHARE::amr
{

/**
 * @brief Stage 1 of the Balsara ADPT divB-free B prolongation.
 *
 * Fills ALL fine faces of a B component from its coarse faces, per component (primal-even
 * direction: exact copy; primal-odd direction: directionalInterp half-point; dual directions:
 * directionalProlongation ±¼ ladder). This is exactly CompositeFieldRefiner run with the magnetic
 * (isMagnetic=true) round-out on — B carries no other special stage-1 gating.
 *
 * Its stage-2 partner is ADPTMagneticRefinePatchStrategy::postprocessRefine: an order-independent
 * divB touch-up that equalizes the 2^d subzone divergences of each coarse zone via a closed-form
 * min-norm correction, making the composite result divB-free exactly regardless of the stage-1
 * order. At order 2 this reproduces the legacy Tóth-Roe operator exactly (derivation §S6b).
 */
template<typename GridLayoutT, typename FieldT, std::size_t order>
using MagneticCompositeRefiner
    = CompositeFieldRefiner<GridLayoutT, FieldT, order, /*isMagnetic=*/true>;


template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>> makeMagneticRefineKernel(int const order)
{
    if (order != 2)
        throw std::runtime_error("makeMagneticRefineKernel: order must be 2 (Linear)");
    return std::make_unique<MagneticCompositeRefiner<GridLayoutT, FieldT, 2>>();
}


} // namespace PHARE::amr


#endif
