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
 * @brief Composable B-field prolongation, restricted to shared (even-normal) faces.
 *
 * B is point-sampled (primal) in its face-normal direction and dual (face-average) in the two
 * tangential directions. The only div-free-compatible, higher-orderable freedom is the tangential
 * dual stencil on the shared faces (copy in normal ⊗ Primitive-B in each tangential dir); the
 * correction is antisymmetric in the child sign ⇒ zero-mean per coarse face ⇒ ∇·B preserved at
 * every order. The interior (odd-normal) fine faces are NOT produced here: they are
 * computed div-free by the Tóth-Roe MagneticRefinePatchStrategy::postprocessRefine, which overwrites
 * them and reads only shared faces + dual children. This is the sharedFacesOnly=true specialization
 * of CompositeFieldRefiner.
 */
template<typename GridLayoutT, typename FieldT, std::size_t order>
using MagneticCompositeRefiner
    = CompositeFieldRefiner<GridLayoutT, FieldT, order, /*sharedFacesOnly=*/true>;


template<typename GridLayoutT, typename FieldT>
std::shared_ptr<IFieldRefineKernel<GridLayoutT, FieldT>>
makeMagneticRefineKernel(int order)
{
    switch (order)
    {
        case 2: return std::make_shared<MagneticCompositeRefiner<GridLayoutT, FieldT, 2>>();
        default: throw std::runtime_error("makeMagneticRefineKernel: order must be 2 (Linear)");
    }
}


} // namespace PHARE::amr


#endif
