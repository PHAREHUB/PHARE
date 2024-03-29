#ifndef PHARE_HDF5_UTILS_HPP
#define PHARE_HDF5_UTILS_HPP

#include "core/utilities/types.hpp"

namespace PHARE::hdf5
{
template<typename T, std::size_t dimension>
inline constexpr auto is_array_dataset
    = (core::is_std_array_v<T, dimension> || core::is_std_array_v<T, 3>);

} // namespace PHARE::hdf5

#endif // PHARE_HDF5_UTILS_HPP
