#ifndef _PHARE_CORE_DATA_FIELD_INITIAZILIZERS_FIELD_USER_INITIALIZER_HPP_
#define _PHARE_CORE_DATA_FIELD_INITIAZILIZERS_FIELD_USER_INITIALIZER_HPP_



#include "core/data/grid/grid_tiles.hpp"

#include <tuple>


namespace PHARE::core
{

class FieldUserFunctionInitializer
{
    template<typename Field, typename GridLayout, typename InitFunction>
    void static _initialize(Field& field, GridLayout const& layout, InitFunction const& init)
    {
        auto const indices = layout.indices(layout.AMRGhostBoxFor(field));
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

public:
    template<typename Field, typename GridLayout, typename InitFunction>
    void static initialize(Field& field, GridLayout const& layout, InitFunction const& init)
    {
        auto constexpr static is_field_tile_set = is_field_tile_set_v<Field>;

        if constexpr (is_field_tile_set)
            for (auto& tile : field())
                _initialize(tile(), tile.layout(), init);
        else
            _initialize(field, layout, init);
    }
};


} // namespace PHARE::core

#endif
