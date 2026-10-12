#ifndef PHARE_PYTHON_PATCH_DATA_HPP
#define PHARE_PYTHON_PATCH_DATA_HPP

#include "pybind_def.hpp"

#include <string>
#include <cstring>
#include <utility>


namespace PHARE::pydata
{
template<typename Data, std::size_t dim>
struct __attribute__((visibility("hidden"))) PatchData
{
    static auto constexpr dimension = dim;
    std::string patchID;
    py_array_t<double> origin{dim};
    py_array_t<int> lower{dim};
    py_array_t<int> upper{dim};
    std::size_t nGhosts;
    Data data;

    PatchData() = default;

    template<typename... Args>
    PatchData(Args&&... args)
        requires std::is_constructible_v<Data, Args&&...>
        : data{std::forward<Args>(args)...}
    {
    }
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
void setPatchDataFromGrid(PatchData& pdata, GridLayout& grid, std::string patchID)
{
    setPatchData(pdata, patchID, *grid.origin(), *grid.AMRBox().lower, *grid.AMRBox().upper);
}

template<typename PatchData, typename Field, typename GridLayout>
void setPatchDataFromField(PatchData& pdata, Field const& field, GridLayout& grid,
                           std::string patchID)
{
    setPatchDataFromGrid(pdata, grid, patchID);
    pdata.nGhosts = GridLayout::options.field_ghost_width;
    pdata.data    = field_as_memory_view(field);
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PATCH_DATA_H*/
