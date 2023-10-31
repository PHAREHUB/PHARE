#ifndef PHARE_TENSORFIELD_HPP
#define PHARE_TENSORFIELD_HPP

#include <cstddef>
#include <string>
#include <array>
#include <vector>
#include <unordered_map>

#include "core/def.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/types.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

namespace PHARE::core
{


template<typename NdArrayImpl, typename PhysicalQuantity, std::size_t rank = 1>
class TensorField
{
private:
    constexpr static std::size_t dimFromRank()
    {
        if constexpr (rank == 1) // Vector field
            return 3;
        else if constexpr (rank == 2) // symmetric 3x3 tensor field
            return 6;
    }

public:
    TensorField()                                     = delete;
    TensorField(TensorField const& source)            = delete;
    TensorField(TensorField&& source)                 = default;
    TensorField& operator=(TensorField const& source) = delete;
    TensorField& operator=(TensorField&& source)      = default;

    using value_type                       = typename NdArrayImpl::type;
    static constexpr std::size_t dimension = NdArrayImpl::dimension;
    using tensor_t                         = typename PhysicalQuantity::template TensorType<rank>;

    TensorField(std::string const& name, tensor_t physQty)
        : name_{name}
        , physQties_{PhysicalQuantity::componentsQuantities(physQty)}
    {
    }
    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    struct TensorFieldProperties
    {
        std::string name;
        typename PhysicalQuantity::Scalar qty;
    };

    using resources_properties = std::array<TensorFieldProperties, dimFromRank()>;

    using field_type = Field<NdArrayImpl, typename PhysicalQuantity::Scalar>;


    NO_DISCARD resources_properties getFieldNamesAndQuantities() const
    {
        return makeResProp_(std::make_index_sequence<dimFromRank()>{});
    }


    void setBuffer(std::string const& bufferName, field_type* field)
    {
        if (auto it = nameToIndex_.find(bufferName); it != std::end(nameToIndex_))
            components_[it->second] = field;
        else
        {
            throw std::runtime_error(
                "TensorField Error - invalid component name, cannot set buffer");
        }
    }




    //! return true if the TensorField can be used to access component data
    NO_DISCARD bool isUsable() const
    {
        return std::all_of(std::begin(components_), std::end(components_),
                           [](auto const& c) { return c != nullptr; });
    }

    NO_DISCARD bool isSettable() const
    {
        return std::all_of(std::begin(components_), std::end(components_),
                           [](auto const& c) { return c == nullptr; });
    }



    void zero()
    {
        if (isUsable())
        {
            for (auto& component : components_)
                component->zero();
        }
        else
        {
            throw std::runtime_error("Error, cannot zero the VecField because it is not usable");
        }
    }


    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD std::string const& name() const { return name_; }



    NO_DISCARD field_type& getComponent(Component component)
    {
        if (isUsable())
        {
            switch (component)
            {
                case Component::X: return *components_[0];
                case Component::Y: return *components_[1];
                case Component::Z: return *components_[2];

                case Component::XX: return *components_[0];
                case Component::XY: return *components_[1];
                case Component::XZ: return *components_[2];
                case Component::YY: return *components_[3];
                case Component::YZ: return *components_[4];
                case Component::ZZ: return *components_[5];
            }
        }
        throw std::runtime_error("Error - TensorField not usable");
    }




    NO_DISCARD field_type const& getComponent(Component component) const
    {
        if (isUsable())
        {
            switch (component)
            {
                case Component::X: return *components_[0];
                case Component::Y: return *components_[1];
                case Component::Z: return *components_[2];

                case Component::XX: return *components_[0];
                case Component::XY: return *components_[1];
                case Component::XZ: return *components_[2];
                case Component::YY: return *components_[3];
                case Component::YZ: return *components_[4];
                case Component::ZZ: return *components_[5];
            }
        }
        throw std::runtime_error("Error - TensorField not usable");
    }



    NO_DISCARD std::string getComponentName(Component component) const
    {
        switch (component)
        {
            case Component::X: return componentNames_[0];
            case Component::Y: return componentNames_[1];
            case Component::Z: return componentNames_[2];
            case Component::XX: return componentNames_[0];
            case Component::XY: return componentNames_[1];
            case Component::XZ: return componentNames_[2];
            case Component::YY: return componentNames_[3];
            case Component::YZ: return componentNames_[4];
            case Component::ZZ: return componentNames_[5];
        }
        throw std::runtime_error("Error - TensorField not usable");
    }

    template<std::size_t... Index>
    NO_DISCARD auto components(std::index_sequence<Index...>) const
    {
        return std::forward_as_tuple((*this)[Index]...); //(*this)[0], (*this)[1], ...);
    }

    template<std::size_t... Index>
    NO_DISCARD auto components() const
    {
        return components(std::make_index_sequence<dimFromRank()>{});
    }
    template<std::size_t... Index>
    NO_DISCARD auto components(std::index_sequence<sizeof...(Index)>)
    {
        return std::forward_as_tuple((*this)[Index]...); //(*this)[0], (*this)[1], ...);
    }

    template<std::size_t... Index>
    NO_DISCARD auto components()
    {
        return components(std::make_index_sequence<dimFromRank()>{});
    }
    NO_DISCARD auto& operator()(Component component) const { return getComponent(component); }
    NO_DISCARD auto& operator()(Component component) { return getComponent(component); }

    NO_DISCARD auto operator()() const { return components(); }
    NO_DISCARD auto operator()() { return components(); }

    NO_DISCARD auto& operator[](std::size_t i) { return *components_[i]; }
    NO_DISCARD auto& operator[](std::size_t i) const { return *components_[i]; }


    void copyData(TensorField const& source)
    {
        if (isUsable() && source.isUsable())
        {
            for (std::size_t i = 0; i < dimFromRank(); ++i)
            {
                components_[i]->copyData(*source.components_[i]);
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




private:
    const std::string name_{"No Name"};

    auto constexpr static N = dimFromRank();
    std::array<typename PhysicalQuantity::Scalar, N> physQties_;

    template<std::size_t... Index>
    resources_properties makeResProp_(std::index_sequence<Index...>) const
    {
        std::array<TensorFieldProperties, sizeof...(Index)> result;
        ((result[Index] = TensorFieldProperties{componentNames_[Index], physQties_[Index]}), ...);
        return result;
    }

    std::array<std::string, N> makeNames_()
    {
        if constexpr (N == 6)
        {
            return {{name_ + "_xx", name_ + "_xy", name_ + "_xz", name_ + "_yy", name_ + "_yz",
                     name_ + "_zz"}};
        }
        else if constexpr (N == 3)
        {
            return {{name_ + "_x", name_ + "_y", name_ + "_z"}};
        }
    }


    template<std::size_t... Index>
    auto makeMap_(std::index_sequence<Index...>) const
    {
        std::unordered_map<std::string, std::size_t> m;
        ((m[componentNames_[Index]] = Index), ...);
        return m;
    }

    const std::array<std::string, N> componentNames_{makeNames_()};
    std::array<field_type*, N> components_{core::ConstArray<field_type*, N>(nullptr)};

    const std::unordered_map<std::string, std::size_t> nameToIndex_{
        makeMap_(std::make_index_sequence<dimFromRank()>{})};
};




template<typename NdArrayImpl, typename PhysicalQuantity>
using SymTensorField = TensorField<NdArrayImpl, PhysicalQuantity, /*rank=*/2>;
} // namespace PHARE::core


#endif /* PHARE_TENSORFIELD_HPP */
