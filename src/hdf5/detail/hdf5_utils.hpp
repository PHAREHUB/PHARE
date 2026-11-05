#ifndef PHARE_HDF5_UTILS_HPP
#define PHARE_HDF5_UTILS_HPP

#include "core/utilities/meta/container_meta.hpp"

namespace PHARE::hdf5
{
template<typename T, std::size_t dimension>
concept IsArrayDataset
    = core::IsStdArray<T> and (std::tuple_size_v<T> == dimension or std::tuple_size_v<T> == 3);

} // namespace PHARE::hdf5

#endif // PHARE_HDF5_UTILS_HPP
