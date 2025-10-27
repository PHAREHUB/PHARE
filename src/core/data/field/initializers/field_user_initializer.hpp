#ifndef _PHARE_CORE_DATA_FIELD_INITIAZILIZERS_FIELD_USER_INITIALIZER_HPP_
#define _PHARE_CORE_DATA_FIELD_INITIAZILIZERS_FIELD_USER_INITIALIZER_HPP_

#include <memory>
#include <tuple>

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"

#include "core/utilities/span.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::core
{
class FieldUserFunctionInitializer
{
public:
    template<typename Field, typename GridLayout>
    void static initialize(Field& field, GridLayout const& layout,
                           initializer::InitFunction<GridLayout::dimension> const& init)
    {
        auto const indices = layout.ghostStartToEndIndices(field, /*includeEnd=*/true);
        auto const coords  = layout.template indexesToCoordVectors</*WithField=*/true>(
            indices, field, [](auto& gridLayout, auto& field_, auto const&... args) {
                return gridLayout.fieldNodeCoordinates(field_, gridLayout.origin(), args...);
            });

        std::shared_ptr<Span<double>> gridPtr // keep grid data alive
            = std::apply([&](auto&... args) { return init(args...); }, coords);
        Span<double>& grid = *gridPtr;

        for (std::size_t cell_idx = 0; cell_idx < indices.size(); cell_idx++)
            std::apply([&](auto&... args) { field(args...) = grid[cell_idx]; }, indices[cell_idx]);

        auto [start_x, end_x] = layout.physicalStartToEnd(field, Direction::X);
        auto [start_y, end_y] = layout.physicalStartToEnd(field, Direction::Y);

        if (field.name() == "mhd_state_B_x")
        {
            end_x -= 1;
            end_y += 1;
            for (std::uint32_t ix = start_x; ix <= end_x; ++ix)
            {
                for (std::uint32_t iy = start_y; iy <= end_y; ++iy)
                {
                    auto index = MeshIndex<GridLayout::dimension>{ix, iy};

                    if (index == MeshIndex<GridLayout::dimension>{0 + 2, 0 + 2}
                        || index == MeshIndex<GridLayout::dimension>{74 + 2, 0 + 2}
                        || index == MeshIndex<GridLayout::dimension>{149 + 2, 50 + 2 + 1}
                        || index == MeshIndex<GridLayout::dimension>{0 + 2, 50 + 2 + 1})
                    {
                        auto coord_x = std::get<0>(coords);
                        std::cout << std::setprecision(16) << "( " << index.str() << ") coord_x--: "
                                  << coord_x(layout.template previous<Direction::Y>(index))
                                  << " coord_x+-: "
                                  << coord_x(layout.template previous<Direction::Y>(
                                         layout.template next<Direction::X>(index)))
                                  << " coord_x-+: " << coord_x(index) << " coord_x++: "
                                  << coord_x(layout.template next<Direction::X>(index)) << "\n";
                    }
                }
            }
        }

        if (field.name() == "mhd_state_B_y")
        {
            end_x -= 1;
            end_y += 1;
            for (std::uint32_t ix = start_x; ix <= end_x; ++ix)
            {
                for (std::uint32_t iy = start_y; iy <= end_y; ++iy)
                {
                    auto index = MeshIndex<GridLayout::dimension>{ix, iy};

                    if (index == MeshIndex<GridLayout::dimension>{0 + 2, 0 + 2}
                        || index == MeshIndex<GridLayout::dimension>{74 + 2 + 1, 0 + 2}
                        || index == MeshIndex<GridLayout::dimension>{149 + 2 + 1, 50 + 2}
                        || index == MeshIndex<GridLayout::dimension>{0 + 2, 50 + 2})
                    {
                        auto coord_y = std::get<1>(coords);
                        std::cout << std::setprecision(16) << "( " << index.str() << ") coord_y--: "
                                  << coord_y(layout.template previous<Direction::X>(index))
                                  << " coord_y+-: "
                                  << coord_y(layout.template previous<Direction::X>(
                                         layout.template next<Direction::Y>(index)))
                                  << " coord_y-+: " << coord_y(index) << " coord_y++: "
                                  << coord_y(layout.template next<Direction::Y>(index)) << "\n";
                    }
                }
            }
        }
    }
};

} // namespace PHARE::core

#endif
