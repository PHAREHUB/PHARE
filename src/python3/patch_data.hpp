#ifndef PHARE_PYTHON_PATCH_DATA_HPP
#define PHARE_PYTHON_PATCH_DATA_HPP


#include "core/utilities/types.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/utilities/point/point.hpp"

#include "amr/resources_manager/amr_utils.hpp"

#include "pybind_def.hpp"


#include <string>
#include <cstring>
#include <utility>
#include <stdexcept>



namespace PHARE::pydata
{
namespace py = pybind11;

template<typename Data, std::size_t dim>
struct __attribute__((visibility("hidden"))) PatchData
{
    static auto constexpr dimension = dim;

    PatchData() = default;

    template<typename... Args>
    PatchData(Args&&... args)
        requires std::is_constructible_v<Data, Args&&...>
        : data{std::forward<Args>(args)...}
    {
    }

    Data data;
    std::string patchID;
    py_array_t<double> origin{dim};
    py_array_t<int> lower{dim};
    py_array_t<int> upper{dim};
    std::size_t nGhosts;
};


template<typename PatchData>
void setPatchData(PatchData& data, std::string const patchID, auto const origin, auto const lower,
                  auto const upper)
{
    std::memcpy(data.lower.request().ptr, lower.data(), PatchData::dimension * sizeof(int));
    std::memcpy(data.upper.request().ptr, upper.data(), PatchData::dimension * sizeof(int));
    std::memcpy(data.origin.request().ptr, origin.data(), PatchData::dimension * sizeof(double));
    data.patchID = patchID;
}

template<typename PatchData, typename GridLayout>
void setPatchDataFromGrid(PatchData& pdata, GridLayout& grid, std::string const& patchID)
{
    setPatchData(pdata, patchID, *grid.origin(), *grid.AMRBox().lower, *grid.AMRBox().upper);
}


template<typename PatchData, typename Field, typename GridLayout>
void setPyPatchDataFromField(PatchData& pdata, Field const& field, GridLayout& grid,
                             std::string patchID)
{
    auto constexpr dimension = PatchData::dimension;
    static_assert(dimension >= 1 and dimension <= 3);

    setPatchDataFromGrid(pdata, grid, patchID);
    pdata.nGhosts = static_cast<std::size_t>(
        GridLayout::nbrGhosts(GridLayout::centering(field.physicalQuantity())[0]));
    pdata.data = field_as_memory_view(field);
}

template<typename Field>
void populatePatchDataFromField(auto& patch_datas, Field const& field, auto const& gridlayout,
                                auto const& patch_id)
{
}

namespace detail
{
    auto static inline accessor = [](auto& el) { return el; };

} // namespace detail

template<typename GridLayout, typename Accessor = decltype(detail::accessor)>
auto getPatchDatasForLevel(auto& hierarchy, auto& model, auto const ilvl, auto& qty,
                           Accessor accessor = detail::accessor)
{
    static constexpr std::size_t dimension = GridLayout::dimension;

    auto& rm = *model.resourcesManager;

    std::vector<PatchData<py_array_t<double>, dimension>> patchDatas;

    auto visit = [&](auto& grid, auto patchID, auto /*ilvl*/) {
        auto&& field = accessor(qty);
        using Field  = std::decay_t<decltype(field)>;

        if constexpr (core::is_field_v<Field>)
            setPyPatchDataFromField(patchDatas.emplace_back(), field, grid, patchID);
        else if constexpr (core::is_field_tile_set_v<Field>)
        {
            std::size_t tidx = 0;
            for (auto const& tile : field())
            {
                auto const tile_id = patchID + "_t" + std::to_string(tidx);
                auto& tile_field   = tile();
                setPyPatchDataFromField(patchDatas.emplace_back(), tile_field, tile.layout(),
                                        tile_id);
                ++tidx;
            }
        }
        else
            static_assert(core::dependent_false_v<GridLayout>);
    };
    amr::visitLevel<GridLayout>(*hierarchy.getPatchLevel(ilvl), rm, visit, qty);

    return patchDatas;
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PATCH_DATA_H*/
