#ifndef PHARE_AMR_RESOURCES_MANAGER_RESOURCES_META_HPP
#define PHARE_AMR_RESOURCES_MANAGER_RESOURCES_META_HPP

#include "core/data/field/field_traits.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"

namespace PHARE::amr
{
/** a field: has a physicalQuantity() */
template<typename ResourceView>
concept IsFieldResource = core::HasPhysicalQuantity<ResourceView>;

/** a tensor field: has components() */
template<typename ResourceView>
concept IsTensorFieldResource = requires(ResourceView view) { view.components(); };

/** particles: can be given a ParticlesPack buffer */
template<typename ResourceView>
concept IsParticlesResource = requires(ResourceView view) {
    view.setBuffer(
        static_cast<core::ParticlesPack<typename ResourceView::particle_array_type>*>(nullptr));
};

template<typename ResourceView>
concept IsResource = IsFieldResource<ResourceView> or IsTensorFieldResource<ResourceView>
                     or IsParticlesResource<ResourceView>;


/** @brief HasRunTimeSubResources: the ResourceView has a runtime list of ResourceViews, like a
 * vector of ResourceViews.
 */
template<typename ResourceView>
concept HasRunTimeSubResources
    = requires(ResourceView view) { view.getRunTimeResourcesViewList(); };


/** @brief HasCompileTimeSubResources: the ResourceView has one or several ResourceViews that
 * can be put in a compile-time list.
 */
template<typename ResourceView>
concept HasCompileTimeSubResources
    = requires(ResourceView view) { view.getCompileTimeResourcesViewList(); };


template<typename ResourceView>
concept HasSubResources
    = HasCompileTimeSubResources<ResourceView> or HasRunTimeSubResources<ResourceView>;

} // namespace PHARE::amr

#endif // PHARE_AMR_RESOURCES_MANAGER_RESOURCES_META_HPP
