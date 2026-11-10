#ifndef PHARE_CORE_DATA_FIELD_ZERO_FIELD_HPP
#define PHARE_CORE_DATA_FIELD_ZERO_FIELD_HPP

#include "core/def.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"

#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>

namespace PHARE::core
{
template<typename FieldT>
class ZeroField
{
public:
    static constexpr std::size_t dimension = FieldT::dimension;

    using type                   = FieldT::type;
    using value_type             = type;
    using physical_quantity_type = FieldT::physical_quantity_type;

    explicit constexpr ZeroField(physical_quantity_type const qty)
        : qty_{qty}
    {
    }

    NO_DISCARD constexpr physical_quantity_type physicalQuantity() const { return qty_; }

    NO_DISCARD std::string name() const { return "zero"; }

    NO_DISCARD constexpr type operator()(auto const&...) const { return type{0}; }

private:
    physical_quantity_type qty_;
};


template<typename VecFieldT>
class ZeroVecField;

template<typename Field_t, typename PhysicalQuantity>
class ZeroVecField<TensorField<Field_t, PhysicalQuantity, 1>>
{
    using vecfield_type = TensorField<Field_t, PhysicalQuantity, 1>;

public:
    static constexpr std::size_t dimension = vecfield_type::dimension;

    using field_type = ZeroField<Field_t>;
    using value_type = vecfield_type::value_type;
    using tensor_t   = vecfield_type::tensor_t;

    explicit constexpr ZeroVecField(tensor_t const qty)
        : qty_{qty}
        , componentQuantities_{PhysicalQuantity::componentsQuantities(qty)}
    {
    }

    NO_DISCARD constexpr field_type operator()(Component const component) const
    {
        return field_type{componentQuantities_[static_cast<std::size_t>(component)]};
    }

    NO_DISCARD constexpr tensor_t physicalQuantity() const { return qty_; }

private:
    tensor_t qty_;
    std::decay_t<decltype(PhysicalQuantity::componentsQuantities(std::declval<tensor_t>()))>
        componentQuantities_;
};


template<typename GridLayout, typename VecFieldT>
NO_DISCARD auto div(GridLayout const&, ZeroVecField<VecFieldT> const&)
{
    return [](MeshIndex<GridLayout::dimension> const) { return 0.; };
}

} // namespace PHARE::core

#endif
