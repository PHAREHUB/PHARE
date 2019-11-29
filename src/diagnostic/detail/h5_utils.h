#ifndef H5_UTILS_H
#define H5_UTILS_H

#include "utilities/types.h"

namespace PHARE::diagnostic::h5
{
template<typename T, std::size_t dimension>
inline constexpr auto is_array_dataset
    = (core::is_std_array_v<T, dimension> || core::is_std_array_v<T, 3>);

} // namespace PHARE::diagnostic::h5

#endif // H5_UTILS_H
