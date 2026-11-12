#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_HPP
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_HPP

#include <array>
#include <utility>
#include <algorithm>
#include <unordered_map>

#include "vecfield_component.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

namespace PHARE
{
namespace core
{
    template<typename Field_t, typename PhysicalQuantity>
    using VecField = TensorField<Field_t, PhysicalQuantity, /*rank=*/1>;



} // namespace core
} // namespace PHARE

#endif
