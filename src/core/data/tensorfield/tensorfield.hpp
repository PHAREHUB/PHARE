#ifndef PHARE_TENSORFIELD_HPP
#define PHARE_TENSORFIELD_HPP

#include <cstddef>
#include <string>
#include <array>
#include <vector>
#include <unordered_map>

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

namespace PHARE::core::detail
{
template<std::size_t rank>
constexpr static std::size_t tensor_field_dim_from_rank()
{
    static_assert(rank > 0 and rank < 3);
    if constexpr (rank == 1) // Vector field
        return 3;
    else if constexpr (rank == 2) // symmetric 3x3 tensor field
        return 6;
}

template<std::size_t rank, typename R = std::array<std::string, tensor_field_dim_from_rank<rank>()>>
R static tensor_field_names(std::string const& name)
{
    if constexpr (rank == 1)
        return {{name + "_x", name + "_y", name + "_z"}};

    else if constexpr (rank == 2)
        return {
            {name + "_xx", name + "_xy", name + "_xz", name + "_yy", name + "_yz", name + "_zz"}};
}

template<typename Field_t, std::size_t N, typename Qtys>
auto static tensor_field_make_fields(std::array<std::string, N> const& names, Qtys const& qtys)
{
    return for_N<N, for_N_R_mode::make_array>([&](auto i) { return Field_t{names[i], qtys[i]}; });
}

} // namespace PHARE::core::detail



namespace PHARE::core
{
template<typename Field_t, typename PhysicalQuantity, std::size_t rank_ = 1>
class TensorField
{
protected:
    auto constexpr static N = detail::tensor_field_dim_from_rank<rank_>();

public:
    static constexpr std::size_t dimension = Field_t::dimension;
    static constexpr std::size_t rank      = rank_;

    using field_type = Field_t;
    using value_type = typename Field_t::type;
    using tensor_t   = typename PhysicalQuantity::template TensorType<rank>;


    TensorField()                                     = delete;
    TensorField(TensorField const& source)            = default;
    TensorField(TensorField&& source)                 = default;
    TensorField& operator=(TensorField const& source) = delete;
    TensorField& operator=(TensorField&& source)      = default;

    TensorField(std::string const& name, tensor_t physQty)
        : qty_{physQty}
        , name_{name}
        , physQties_{PhysicalQuantity::componentsQuantities(physQty)}
        , componentNames_{detail::tensor_field_names<rank>(name)}
        , components_{detail::tensor_field_make_fields<Field_t>(componentNames_, physQties_)}
        , nameToIndex_{makeMap_(std::make_index_sequence<N>{})}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    void setBuffer(std::nullptr_t ptr)
    {
        for_N<N, for_N_R_mode::forward_tuple>([&](auto i) { components_[i].setBuffer(nullptr); });
    }

    template<typename Fields>
    void setBuffer(Fields* const fields)
    {
        if (!fields)
            throw std::runtime_error("use other fn");
        for_N<N>([&](auto i) { components_[i].setBuffer(&(*fields)[i]); });
    }


    //! return true if the TensorField can be used to access component data
    NO_DISCARD bool isUsable() const
    {
        return std::all_of(std::begin(components_), std::end(components_),
                           [](auto const& c) { return c.isUsable(); });
    }

    NO_DISCARD bool isSettable() const
    {
        return std::all_of(std::begin(components_), std::end(components_),
                           [](auto const& c) { return c.isSettable(); });
    }



    void zero()
    {
        if (!isUsable())
            throw std::runtime_error("Error, cannot zero the TensorField because it is not usable");

        for (auto& component : components_)
            component.zero();
    }


    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD std::string const& name() const { return name_; }




    void _check() const
    {
        if (!isUsable())
            throw std::runtime_error("Error - TensorField not usable");
    }

    NO_DISCARD field_type& getComponent(Component component)
    {
        _check();
        return components_[_get_index_for(component)];
    }

    NO_DISCARD field_type const& getComponent(Component component) const
    {
        _check();
        return components_[_get_index_for(component)];
    }


    NO_DISCARD std::string getComponentName(Component component) const
    {
        return componentNames_[_get_index_for(component)];
    }


    template<std::size_t... Index>
    NO_DISCARD auto components(std::index_sequence<Index...>) const
    {
        return std::forward_as_tuple((*this)[Index]...);
    }
    NO_DISCARD auto components() const { return components(std::make_index_sequence<N>{}); }


    template<std::size_t... Index>
    NO_DISCARD auto components(std::index_sequence<Index...>)
    {
        return std::forward_as_tuple((*this)[Index]...);
    }
    NO_DISCARD auto components() { return components(std::make_index_sequence<N>{}); }


    NO_DISCARD auto& operator()(Component component) const { return getComponent(component); }
    NO_DISCARD auto& operator()(Component component) { return getComponent(component); }


    NO_DISCARD auto operator()() const { return components(); }
    NO_DISCARD auto operator()() { return components(); }


    NO_DISCARD auto& operator[](std::size_t i) { return components_[i]; }
    NO_DISCARD auto& operator[](std::size_t i) const { return components_[i]; }


    void copyData(TensorField const& source)
    {
        if (isUsable() && source.isUsable())
        {
            for (std::size_t i = 0; i < N; ++i)
            {
                components_[i].copyData(source.components_[i]);
            }
        }
        else
        {
            throw std::runtime_error("Error, unusable TensorField, cannot copyData");
        }
    }


    NO_DISCARD auto begin() { return std::begin(components_); }
    NO_DISCARD auto cbegin() const { return std::cbegin(components_); }
    NO_DISCARD auto end() { return std::end(components_); }
    NO_DISCARD auto cend() const { return std::cend(components_); }

    NO_DISCARD auto& componentNames() const { return componentNames_; }
    NO_DISCARD auto& physicalQuantity() const { return qty_; }
    NO_DISCARD auto constexpr static size() { return N; }

private:
    auto static _get_index_for(Component component)
    {
        auto val = static_cast<std::underlying_type_t<Component>>(component);
        if constexpr (rank == 1)
            return val;
        else if constexpr (rank == 2)
            return val - detail::tensor_field_dim_from_rank<1>();
    }


    template<std::size_t... Index>
    auto makeMap_(std::index_sequence<Index...>) const
    {
        std::unordered_map<std::string, std::size_t> m;
        ((m[componentNames_[Index]] = Index), ...);
        return m;
    }



    tensor_t qty_;
    std::string const name_{"No Name"};
    std::array<typename PhysicalQuantity::Scalar, N> physQties_;
    std::array<std::string, N> const componentNames_;
    std::array<field_type, N> components_;
    std::unordered_map<std::string, std::size_t> const nameToIndex_;
};



template<typename Field_t, typename PhysicalQuantity>
using SymTensorField = TensorField<Field_t, PhysicalQuantity, /*rank=*/2>;


} // namespace PHARE::core


#endif /* PHARE_TENSORFIELD_HPP */
