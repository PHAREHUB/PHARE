#ifndef PHARE_REFINEMENT_CONFIG_HPP
#define PHARE_REFINEMENT_CONFIG_HPP

#include "initializer/data_provider.hpp"
#include "amr/data/field/refine/field_refiner_kernel.hpp"

#include <stdexcept>
#include <string>

namespace PHARE::amr
{

/**
 * @brief Runtime selection of the field-refinement order.
 *
 * The refine-op members of the messengers are built from the runtime kernels
 * (makeRefineKernel / makeMagneticRefineKernel), whose stencil is selected by this order.
 * FieldRefinementOrder::Linear (order 2) is the only supported value for now; this struct is the
 * extension point for higher orders. Distinct from the EXISTING particle split-operator template
 * param named RefinementParams (MessengerFactory / HybridHybridMessengerStrategy) — do not
 * conflate.
 */
struct RefinementConfig
{
    FieldRefinementOrder order = FieldRefinementOrder::Linear;

    //! Read the optional field-refinement selection from the dict. Absent "AMR"/"refinement"/
    //! "order" nodes at any level of the path ⇒ the RefinementConfig default order. Only order 2
    //! (Linear) is supported; this is the one place a raw dict value is validated against it.
    RefinementConfig static FROM(PHARE::initializer::PHAREDict const& dict)
    {
        PHARE::amr::RefinementConfig config;
        auto const rawOrder = cppdict::get_value(dict, "simulation/AMR/refinement/order", int{2});
        if (rawOrder != 2)
            throw std::runtime_error("unsupported field refinement order: "
                                     + std::to_string(rawOrder) + " (only 2 is supported)");
        config.order = static_cast<FieldRefinementOrder>(rawOrder);
        return config;
    }
};

} // namespace PHARE::amr

#endif
