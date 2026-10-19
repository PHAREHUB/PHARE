#ifndef DEFAULT_TAGGER_STRATEGY_H
#define DEFAULT_TAGGER_STRATEGY_H

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "initializer/data_provider.hpp"

#include "tagger_strategy.hpp"
#include "tagging_quantities.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>
#include <algorithm>

namespace PHARE::amr
{
template<typename Model>
class DefaultTaggerStrategy : public TaggerStrategy<Model>
{
    using gridlayout_type           = Model::gridlayout_type;
    static auto constexpr dimension = Model::dimension;
    using components_type           = std::vector<tagging::tagged_component_t<Model>>;


public:
    DefaultTaggerStrategy(initializer::PHAREDict const& dict)
        : quantities_{tagging::parseQuantities(dict)}
    {
    }
    void tag(Model& model, gridlayout_type const& layout, int* tags) const override;

private:
    template<typename Diff>
    static double maxOverComponents_(components_type const& components, Diff&& diff)
    {
        auto values = diff(components[0]);
        auto crit   = values[0];
        for (std::size_t iComp = 0; iComp < components.size(); ++iComp)
        {
            if (iComp > 0)
                values = diff(components[iComp]);
            for (auto const value : values)
                if (crit < value)
                    crit = value;
        }
        return crit;
    }

    std::vector<tagging::QuantityThreshold> quantities_;
};

template<typename Model>
void DefaultTaggerStrategy<Model>::tag(Model& model, gridlayout_type const& layout, int* tags) const
{
    auto const quantities = tagging::resolveQuantities(model, quantities_);

    auto const exceeds = [&](auto const& criterion) {
        return std::any_of(std::begin(quantities), std::end(quantities), [&](auto const& quantity) {
            return criterion(quantity.first) > quantity.second;
        });
    };

    // we loop on cell indexes for all qties regardless of their centering
    auto const& start_x = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::X);

    // override end_x because the tag buffer does not have ghost cells
    // and physicalEnd will account for ghost cells
    auto const& end_x = layout.nbrCells()[0] - 1;

    // SAMRAI tags int* buffer is FORTRAN ordering so we set false to the view
    bool constexpr c_ordering = false;
    auto tagsv = core::NdArrayView<dimension, int, c_ordering>(tags, layout.nbrCells());

    if constexpr (dimension == 1)
    {
        // at interporder 1 we choose not to tag the last patch cell since
        // the 5 points stencil may go beyond the last ghost node.
        // for interp order 2 and 3 this is ok
        auto constexpr doLastCell = gridlayout_type::options.field_ghost_width > 2;
        std::size_t oneOrZero     = doLastCell ? 1 : 0;
        for (auto iCell = 0u, ix = start_x; iCell < end_x + oneOrZero; ++ix, ++iCell)
        {
            auto const criterion = [&](components_type const& components) {
                double sum = 0.;
                for (auto const& F : components)
                {
                    auto Favg   = 0.2 * (F(ix - 2) + F(ix - 1) + F(ix) + F(ix + 1) + F(ix + 2));
                    auto Favgp1 = 0.2 * (F(ix - 1) + F(ix) + F(ix + 1) + F(ix + 2) + F(ix + 3));
                    auto criter = std::abs(Favgp1 - Favg) / (1 + std::abs(Favg));
                    sum += criter * criter;
                }
                return std::sqrt(sum);
            };

            tagsv(iCell) = exceeds(criterion) ? 1 : 0;
        }
    }
    if constexpr (dimension == 2)
    {
        auto const& start_y
            = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::Y);

        auto const& end_y = layout.nbrCells()[1] - 1;

        for (auto iTag_x = 0u, ix = start_x; iTag_x <= end_x; ++ix, ++iTag_x)
        {
            for (auto iTag_y = 0u, iy = start_y; iTag_y <= end_y; ++iy, ++iTag_y)
            {
                auto field_diff = [&](auto const& F) //
                {
                    auto const delta_2x = std::abs(F(ix + 2, iy) - F(ix, iy));
                    auto const delta_2y = std::abs(F(ix, iy + 2) - F(ix, iy));
                    auto const delta_x  = std::abs(F(ix + 1, iy) - F(ix, iy));
                    auto const delta_y  = std::abs(F(ix, iy + 1) - F(ix, iy));

                    auto const criter_x = delta_2x / (1 + delta_x);
                    auto const criter_y = delta_2y / (1 + delta_y);

                    return std::array{criter_x, criter_y};
                };

                auto const criterion = [&](components_type const& components) {
                    return maxOverComponents_(components, field_diff);
                };

                tagsv(iTag_x, iTag_y) = exceeds(criterion) ? 1 : 0;
            }
        }
    }
    if constexpr (dimension == 3)
    {
        auto const& start_y
            = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::Y);
        auto const& start_z
            = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::Z);

        auto const& end_y = layout.nbrCells()[1] - 1;
        auto const& end_z = layout.nbrCells()[2] - 1;

        for (auto iTag_x = 0u, ix = start_x; iTag_x <= end_x; ++ix, ++iTag_x)
        {
            for (auto iTag_y = 0u, iy = start_y; iTag_y <= end_y; ++iy, ++iTag_y)
            {
                for (auto iTag_z = 0u, iz = start_z; iTag_z <= end_z; ++iz, ++iTag_z)
                {
                    auto const field_diff = [&](auto const& F) {
                        auto const delta_2x = std::abs(F(ix + 2, iy, iz) - F(ix, iy, iz));
                        auto const delta_2y = std::abs(F(ix, iy + 2, iz) - F(ix, iy, iz));
                        auto const delta_2z = std::abs(F(ix, iy, iz + 2) - F(ix, iy, iz));
                        auto const delta_x  = std::abs(F(ix + 1, iy, iz) - F(ix, iy, iz));
                        auto const delta_y  = std::abs(F(ix, iy + 1, iz) - F(ix, iy, iz));
                        auto const delta_z  = std::abs(F(ix, iy, iz + 1) - F(ix, iy, iz));

                        auto const criter_x = delta_2x / (1 + delta_x);
                        auto const criter_y = delta_2y / (1 + delta_y);
                        auto const criter_z = delta_2z / (1 + delta_z);

                        return std::array{criter_x, criter_y, criter_z};
                    };

                    auto const criterion = [&](components_type const& components) {
                        return maxOverComponents_(components, field_diff);
                    };

                    tagsv(iTag_x, iTag_y, iTag_z) = exceeds(criterion) ? 1 : 0;
                }
            }
        }
    }
}
} // namespace PHARE::amr

#endif // DEFAULT_HYBRID_TAGGER_STRATEGY_H
