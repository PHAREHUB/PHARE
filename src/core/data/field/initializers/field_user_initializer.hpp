#ifndef _PHARE_CORE_DATA_FIELD_INITIAZILIZERS_FIELD_USER_INITIALIZER_HPP_
#define _PHARE_CORE_DATA_FIELD_INITIAZILIZERS_FIELD_USER_INITIALIZER_HPP_

#include "core/utilities/span.hpp"
#include "initializer/data_provider.hpp"
#include "core/utilities/point/point.hpp"

#include <tuple>
#include <memory>
#include <cassert>

namespace PHARE::core
{
class FieldUserFunctionInitializer
{
public:
    template<typename Field, typename GridLayout>
    void static initialize(Field& field, GridLayout const& layout,
                           initializer::InitFunction<GridLayout::dimension> const& init)
    {
        auto const indices = layout.indicis(layout.AMRGhostBoxFor(field));
        auto const coords  = layout.template indexesToCoordVectors</*WithField=*/true>(
            indices, field, [](auto& gridLayout, auto& field_, auto const&... args) {
                return gridLayout.fieldNodeCoordinates(field_, args...);
            });

        std::shared_ptr<Span<double>> gridPtr // keep grid data alive
            = std::apply([&](auto&... args) { return init(args...); }, coords);
        Span<double>& grid = *gridPtr;

        for (std::size_t cell_idx = 0; cell_idx < indices.size(); cell_idx++)
            std::apply(
                [&](auto&... args) { field(layout.AMRToLocal(Point{args...})) = grid[cell_idx]; },
                indices[cell_idx]);
    }
};

} // namespace PHARE::core

#endif
