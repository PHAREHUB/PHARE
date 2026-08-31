#ifndef PHARE_MAGNETIC_COMPOSITE_REFINER_HPP
#define PHARE_MAGNETIC_COMPOSITE_REFINER_HPP


#include "phare_mpi.hpp" // IWYU pragma: keep

#include "field_refiner_kernel.hpp"
#include "composite_field_refiner.hpp"

#include <memory>


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
 * order. At order 2 this reproduces the legacy Tóth-Roe operator exactly.
 */
template<typename GridLayoutT, typename FieldT, std::size_t order>
using MagneticCompositeRefiner
    = CompositeFieldRefiner<GridLayoutT, FieldT, order, /*isMagnetic=*/true>;


// Single-enumerator dispatch, same as makeRefineKernel above: the order is validated once in
// RefinementConfig::FROM and the enum type carries it from there.
template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>>
makeMagneticRefineKernel([[maybe_unused]] FieldRefinementOrder const order)
{
    return std::make_unique<MagneticCompositeRefiner<
        GridLayoutT, FieldT, static_cast<std::size_t>(FieldRefinementOrder::Linear)>>();
}


} // namespace PHARE::amr


#endif
