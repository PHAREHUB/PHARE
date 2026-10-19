#ifndef PHARE_CORE_NUMERICS_EXTERNAL_FIELD_PROJECTION_HPP
#define PHARE_CORE_NUMERICS_EXTERNAL_FIELD_PROJECTION_HPP

#include "core/def.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"

namespace PHARE::core
{
template<typename GridLayout, auto src, auto dst, typename VecFieldT>
NO_DISCARD auto projectB0(VecFieldT const* B0, Component const component,
                          MeshIndex<GridLayout::dimension> const& index)
{
    using value_type = VecFieldT::value_type;

    if (B0 == nullptr)
        return value_type{0};

    return GridLayout::template project<src, dst>((*B0)(component), index);
}

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_EXTERNAL_FIELD_PROJECTION_HPP
