#ifndef PHARE_COMPONENT_WISE_TAGGER_STRATEGY_HPP
#define PHARE_COMPONENT_WISE_TAGGER_STRATEGY_HPP

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
#include <utility>
#include <vector>

namespace PHARE::amr
{
template<typename Derived, typename Model>
class ComponentWiseTaggerStrategy : public TaggerStrategy<Model>
{
public:
    using gridlayout_type           = Model::gridlayout_type;
    static constexpr auto dimension = Model::dimension;
    using recentered_field_type
        = core::RecenteredField<gridlayout_type, tagging::tagged_field_t<Model>>;

    explicit ComponentWiseTaggerStrategy(initializer::PHAREDict const& dict)
        : quantities_{tagging::parseQuantities(dict)}
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

    void tag(Model& model, gridlayout_type const& layout, int* tags) const final;

private:
    static constexpr std::uint32_t ghostWidth_()
    {
        return static_cast<std::uint32_t>(gridlayout_type::options.field_ghost_width);
    }

    static constexpr std::uint32_t margin_()
    {
        auto constexpr reach  = Derived::stencil_reach;
        auto constexpr ghosts = ghostWidth_();
        return reach > ghosts ? reach - ghosts : 0;
    }

    Derived const& derived_() const { return static_cast<Derived const&>(*this); }

    std::vector<tagging::QuantityThreshold> quantities_;
};


template<typename Derived, typename Model>
void ComponentWiseTaggerStrategy<Derived, Model>::tag(Model& model, gridlayout_type const& layout,
                                                      int* tags) const
{
    auto const allDual = core::ConstArray<core::QtyCentering, dimension>(core::QtyCentering::dual);

    std::vector<std::pair<std::vector<recentered_field_type>, double>> quantities;
    for (auto const& [components, threshold] : tagging::resolveQuantities(model, quantities_))
    {
        auto& [fields, quantityThreshold] = quantities.emplace_back();
        quantityThreshold                 = threshold;
        for (auto const& component : components)
            fields.emplace_back(component.get(), allDual);
    }

    constexpr bool hasCombine
        = requires(std::span<double const> indicators) { Derived::combine(indicators); };
    std::array<double, tagging::max_tagged_components_v<Model>> indicators;

    auto const nbrCells = layout.nbrCells();
    std::array<std::uint32_t, dimension> start;
    for (std::size_t d = 0; d < dimension; ++d)
        start[d]
            = layout.physicalStartIndex(core::QtyCentering::dual, static_cast<core::Direction>(d));

    auto constexpr margin = margin_();
    auto const nbrTags    = core::product(nbrCells);

    for (std::size_t i = 0; i < nbrTags; ++i)
    {
        std::array<std::uint32_t, dimension> cell;
        bool inBorder = false;
        auto rest     = i;
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto const t = static_cast<std::uint32_t>(rest % nbrCells[d]);
            rest /= nbrCells[d];
            if constexpr (margin > 0)
                inBorder = inBorder or t < margin or t + margin >= nbrCells[d];
            cell[d] = t + start[d];
        }

        auto const passes = [&](auto const& quantity) {
            auto const& [fields, threshold] = quantity;
            if constexpr (hasCombine)
            {
                for (std::size_t c = 0; c < fields.size(); ++c)
                    indicators[c] = derived_().indicator(fields[c], cell);
                return Derived::combine(std::span<double const>(indicators.data(), fields.size()))
                       > threshold;
            }
            else
                return core::any(fields, [&](auto const& field) {
                    return derived_().indicator(field, cell) > threshold;
                });
        };
        tags[i] = !inBorder and core::any(quantities, passes) ? 1 : 0;
    }
}

} // namespace PHARE::amr

#endif
