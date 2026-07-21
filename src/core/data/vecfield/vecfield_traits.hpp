#ifndef PHARE_CORE_DATA_VECFIELD_TRAITS_HPP
#define PHARE_CORE_DATA_VECFIELD_TRAITS_HPP

#include "core/data/tensorfield/tensorfield_traits.hpp"

namespace PHARE::core
{
template<typename T>
concept IsVecField = PHARE::core::IsTensorField<T> || (T::N == 3);
}
#endif // PHARE_CORE_DATA_VECFIELD_TRAITS_HPP
