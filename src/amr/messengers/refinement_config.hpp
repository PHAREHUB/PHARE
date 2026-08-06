#ifndef PHARE_REFINEMENT_CONFIG_HPP
#define PHARE_REFINEMENT_CONFIG_HPP

#include "initializer/data_provider.hpp"

#include <stdexcept>
#include <string>

namespace PHARE::amr
{

/**
 * @brief Runtime selection of the field-refinement order.
 *
 * The refine-op members of the messengers are built from the runtime kernels
 * (makeRefineKernel / makeMagneticRefineKernel), whose stencil is selected by this order.
 * order == 2 (Linear) is the only supported value for now; this struct is the extension point
 * for higher orders. Distinct from the EXISTING particle split-operator template param named
 * RefinementParams (MessengerFactory / HybridHybridMessengerStrategy) — do not conflate.
 */
struct RefinementConfig
{
    int order = 2;

    //! Read the optional field-refinement selection from the dict. Absent keys ⇒ the
    //! RefinementConfig default order. Only order 2 (Linear) is supported.
    RefinementConfig static FROM(PHARE::initializer::PHAREDict const& dict)
    {
        PHARE::amr::RefinementConfig config;
        auto const& refinement = dict["simulation"]["AMR"]["refinement"];
        if (refinement.contains("order"))
            config.order = refinement["order"].template to<int>();
        if (config.order != 2)
            throw std::runtime_error("unsupported field refinement order: "
                                     + std::to_string(config.order) + " (only 2 is supported)");
        return config;
    }
};

} // namespace PHARE::amr

#endif
