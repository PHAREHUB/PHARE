#ifndef PHARE_PYTHON_PATCH_DATA_H
#define PHARE_PYTHON_PATCH_DATA_H

#include "pybind_def.h"

namespace PHARE::pydata
{
template<typename Data, std::size_t dim>
struct PatchData
{
    static auto constexpr dimension = dim;
    std::string patchID;
    std::string origin;
    py_array_t<std::size_t> lower{dim};
    py_array_t<std::size_t> upper{dim};
    std::size_t nGhosts;
    Data data;

    PatchData() = default;

    template<typename... Args>
    PatchData(Args&&... args)
        : data{std::forward<Args...>(args...)}
    {
    }
};


template<typename PatchData>
void setPatchData(PatchData& data, std::string patchID, std::string origin,
                  std::array<std::size_t, PatchData::dimension> lower,
                  std::array<std::size_t, PatchData::dimension> upper)
{
    constexpr std::size_t bytes = PatchData::dimension * sizeof(size_t);
    std::memcpy(data.lower.request().ptr, lower.data(), bytes);
    std::memcpy(data.upper.request().ptr, upper.data(), bytes);
    data.patchID = patchID;
    data.origin  = origin;
}

template<typename PatchData, typename GridLayout>
void setPatchDataFromGrid(PatchData& pdata, GridLayout& grid, std::string patchID)
{
    setPatchData(pdata, patchID, grid.origin().str(),
                 grid.AMRBox().lower.template toArray<std::size_t>(),
                 grid.AMRBox().upper.template toArray<std::size_t>());
}

template<typename PatchData, typename Field, typename GridLayout>
void setPatchDataFromField(PatchData& pdata, Field const& field, GridLayout& grid,
                           std::string patchID)
{
    setPatchDataFromGrid(pdata, grid, patchID);
    pdata.nGhosts = static_cast<std::size_t>(
        GridLayout::nbrGhosts(GridLayout::centering(field.physicalQuantity())[0]));
    pdata.data.assign(field.data(), field.data() + field.size());
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PATCH_DATA_H*/
