#ifndef PHARE_PER_QUANTITY_TAGGER_STRATEGY_HPP
#define PHARE_PER_QUANTITY_TAGGER_STRATEGY_HPP

#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/field/recentered_field.hpp"

#include "initializer/data_provider.hpp"

#include "tagger_strategy.hpp"
#include "tagging_quantities.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace PHARE::amr
{
/**
 * @brief CRTP base for tagging based on a single criterion applied to a list of quantities.
 *
 * The criterion is made of an indicator function (the same for all quantities), and a per-quantity
 * threshold. The criterion is met in a cell if at least one quantity has its indicator exceeding
 * its threshold.
 *
 * This class expects its child to implement the indicator with the following interface:
 *
 * @code
 * double indicator(std::span<recentered_field_type const> components,
 *                  std::array<std::uint32_t, dimension> const& cell,
 *                  gridlayout_type const& layout) const;
 * @endcode
 *
 * where `components` are the components of one quantity, in component order (a single one for a
 * scalar quantity), each viewed at cell centres, `cell` the index of the cell in the field arrays,
 * the first ghost cell having index 0, and `layout` the layout of the patch being tagged. The
 * indicator returns a single value for the quantity, compared to the quantity's threshold.
 *
 * The CRTP child is also expected to provide a `static constexpr std::uint32_t stencil_reach`, the
 * largest number of cells the indicator reads away from `cell`, to avoid out-of-bounds access in
 * case the ghost layer width is not large enough. In this case, a margin at the patch border is
 * left with untagged cells.
 *
 * @tparam Derived The CRTP child class
 * @tparam Model The physical model
 */
template<typename Derived, typename Model>
class PerQuantityTaggerStrategy : public TaggerStrategy<Model>
{
public:
    using gridlayout_type       = Model::gridlayout_type;
    using recentered_field_type = core::RecenteredField<tagging::tagged_field_t<Model>>;

    static constexpr auto dimension = Model::dimension;

    /**
     * @brief Constructor
     *
     * @param dict config dict at path "simulation/AMR/refinement/tagging"
     */
    explicit PerQuantityTaggerStrategy(initializer::PHAREDict const& dict)
        : quantities_{tagging::resolveQuantities<Model>(dict)}
    {
        auto constexpr reach  = Derived::stencil_reach;
        auto constexpr ghosts = ghostWidth_();
        if constexpr (reach > ghosts)
        {
            PHARE_LOG_LINE_SS("tagging: criterion stencil reach "
                              << reach << " exceeds the field ghost width " << ghosts
                              << ", the outer " << reach - ghosts
                              << " cell layers of each patch are not tagged");
        }
    }

    /**
     * @brief The main tag function, called once per patch by ConcreteTagger::tag()
     * (amr/tagging/concrete_tagger.hpp), after the model has been set on that patch.
     *
     * @param model The physical model, set on the patch.
     * @param layout The grid layout of the current patch.
     * @param tags The array of tag values, one per cell of the patch, in Fortran ordering.
     */
    void tag(Model& model, gridlayout_type const& layout, int* tags) const final;

private:
    /// Bundles the components' recentered views of a single quantity, along with its threshold.
    struct TaggedQuantityRecenteredComponents
    {
        std::vector<recentered_field_type>
            recenteredFields; ///< recentered views of the components' fields.
        double threshold;     ///< the criterion's threshold for this quantity.
    };

    /// Retrieve the ghost layer width
    static constexpr std::uint32_t ghostWidth_()
    {
        return static_cast<std::uint32_t>(gridlayout_type::options.field_ghost_width);
    }

    /// Return the layer size of untaggable physical cells near the patch border, because of the
    /// stencil size of the indicator function.
    static constexpr std::uint32_t margin_()
    {
        auto constexpr reach  = Derived::stencil_reach;
        auto constexpr ghosts = ghostWidth_();
        return reach > ghosts ? reach - ghosts : 0;
    }

    /// CRTP downcast helper.
    Derived const& derived_() const { return static_cast<Derived const&>(*this); }

    /// The list of each quantity's model accessor (also contains the threshold)
    std::vector<tagging::TaggedQuantityAccessor<Model>> quantities_;
};


template<typename Derived, typename Model>
void PerQuantityTaggerStrategy<Derived, Model>::tag(Model& model, gridlayout_type const& layout,
                                                    int* tags) const
{
    auto const allDual = core::ConstArray<core::QtyCentering, dimension>(core::QtyCentering::dual);

    // build the recentered views of all quantities' components
    std::vector<TaggedQuantityRecenteredComponents> quantities;
    for (auto const& [components, threshold] : quantities_)
    {
        auto& quantity     = quantities.emplace_back();
        quantity.threshold = threshold;
        for (auto const& component : components(model))
            quantity.recenteredFields.emplace_back(component, allDual, layout);
    }

    auto const nbrCells = layout.nbrCells();
    std::array<std::uint32_t, dimension> start;
    for (std::size_t d = 0; d < dimension; ++d)
        start[d]
            = layout.physicalStartIndex(core::QtyCentering::dual, static_cast<core::Direction>(d));

    auto constexpr margin = margin_();
    auto const nbrTags    = core::product(nbrCells);

    for (std::size_t tagId = 0; tagId < nbrTags; ++tagId)
    {
        // determine at the same time the indices of the cell corresponding to tagId, and if the
        // cell is inside the untaggable margin. tagId is a flat index in Fortran ordering (first
        // index fastest), decoded here one direction at a time.
        std::array<std::uint32_t, dimension> cellId;
        bool inUntaggableMargin = false;
        auto rest               = tagId;
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto const t = static_cast<std::uint32_t>(rest % nbrCells[d]);
            rest /= nbrCells[d];
            if constexpr (margin > 0)
                inUntaggableMargin = inUntaggableMargin or t < margin or t + margin >= nbrCells[d];
            cellId[d] = t + start[d];
        }

        // lambda expressing the tagging criterion for a single tagged quantity
        auto const passes = [&](auto const& quantity) {
            auto const& [recenteredFields, threshold] = quantity;
            return derived_().indicator(std::span<recentered_field_type const>(recenteredFields),
                                        cellId, layout)
                   > threshold;
        };

        // this sets directly a null tag if the cell is in the untaggable margin, without evaluating
        // the lambda above, so there is no out-of-bounds memory access.
        tags[tagId] = !inUntaggableMargin and core::any(quantities, passes) ? 1 : 0;
    }
}

} // namespace PHARE::amr

#endif
