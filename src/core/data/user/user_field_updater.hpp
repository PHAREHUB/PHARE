#ifndef PHARE_CORE_DATA_USER_USER_FIELD_UPDATER_HPP
#define PHARE_CORE_DATA_USER_USER_FIELD_UPDATER_HPP

#include "core/utilities/span.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/space_time_function.hpp"

#include <array>
#include <tuple>
#include <memory>
#include <cassert>

namespace PHARE::core
{
/**
 * @brief Fills a field from a user function of space and time.
 *
 */
class UserFieldUpdater
{
public:
    template<typename Field, typename GridLayout>
    void static update(Field& field, GridLayout const& layout,
                       SpaceTimeFunction<GridLayout::dimension> const& f, double time)
    {
        auto const indices = layout.indices(layout.AMRGhostBoxFor(field));
        auto const coords  = layout.template indexesToCoordVectors</*WithField=*/true>(
            indices, field, [](auto& gridLayout, auto& field_, auto const&... args) {
                return gridLayout.fieldNodeCoordinates(field_, args...);
            });

        std::shared_ptr<Span<double>> gridPtr // keep grid data alive
            = std::apply([&](auto const&... xyz) { return f(make_span(xyz)..., time); }, coords);
        Span<double>& grid = *gridPtr;

        // a user function returning the wrong number of values would be read out of bounds
        assert(grid.size() == indices.size());

        for (std::size_t cell_idx = 0; cell_idx < indices.size(); cell_idx++)
            std::apply(
                [&](auto&... args) { field(layout.AMRToLocal(Point{args...})) = grid[cell_idx]; },
                indices[cell_idx]);
    }

    template<typename VecField, typename GridLayout>
    void static update(
        VecField& vecfield, GridLayout const& layout,
        std::array<SpaceTimeFunction<GridLayout::dimension>, VecField::size()> const& funcs,
        double time)
    {
        for (std::size_t i = 0; i < VecField::size(); ++i)
            update(vecfield[i], layout, funcs[i], time);
    }
};
} // namespace PHARE::core

#endif // PHARE_CORE_DATA_USER_USER_FIELD_UPDATER_HPP
