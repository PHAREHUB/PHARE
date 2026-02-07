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


    template<typename VecField, typename = tryToInstanciate<typename VecField::field_type>>
    void average(VecField const& vf1, VecField const& vf2, VecField& Vavg)
    {
        average(vf1.getComponent(Component::X), vf2.getComponent(Component::X),
                Vavg.getComponent(Component::X));

        average(vf1.getComponent(Component::Y), vf2.getComponent(Component::Y),
                Vavg.getComponent(Component::Y));

        average(vf1.getComponent(Component::Z), vf2.getComponent(Component::Z),
                Vavg.getComponent(Component::Z));
    }

} // namespace core
} // namespace PHARE

#endif
