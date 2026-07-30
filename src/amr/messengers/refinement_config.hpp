#ifndef PHARE_REFINEMENT_CONFIG_HPP
#define PHARE_REFINEMENT_CONFIG_HPP

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
};

} // namespace PHARE::amr

#endif
