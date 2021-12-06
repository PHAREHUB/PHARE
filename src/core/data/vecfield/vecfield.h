#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_H
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_H

#include <array>
#include <utility>
#include <algorithm>
#include <unordered_map>

#include "core/def.h"
#include "core/utilities/types.h"
#include "core/utilities/meta/meta_utilities.h"

#include "vecfield_component.h"
#include "core/data/field/field.h"

namespace PHARE::core
{
template<typename Field>
class VecFieldView
{
public:
    using This                             = VecFieldView<Field>;
    using field_type                       = Field;
    using Components                       = std::array<Field, 3>;
    static constexpr std::size_t dimension = Field::dimension;

    template<typename VecField>
    static constexpr Components make_components_from(VecField& vecfield)
    {
        return core::generate(
            [&](auto& field_) {
                auto& field = core::deref(field_); // might be a pointer
                return field_type{field.data(), field.shape(), field.physicalQuantity()};
            },
            vecfield.components_);
    }


    VecFieldView(Components components)
        : components_{components}
    {
    }


    template<typename VecField>
    VecFieldView(VecField& vecfield)
        : components_{make_components_from(vecfield)}
    {
    }

    VecFieldView()                    = default;
    VecFieldView(VecFieldView const&) = default;
    VecFieldView(VecFieldView&&)      = default;

    auto& getComponent(Component component) const _PHARE_ALL_FN_
    {
        switch (component)
        {
            case Component::X: return core::deref(components_[0]);
            case Component::Y: return core::deref(components_[1]);
            case Component::Z: return core::deref(components_[2]);
        }
        throw_runtime_error("Error - VecFieldView not usable");
        return core::deref(components_[0]); // hax
    }


    auto& getComponent(Component component) _PHARE_ALL_FN_
    {
        switch (component)
        {
            case Component::X: return core::deref(components_[0]);
            case Component::Y: return core::deref(components_[1]);
            case Component::Z: return core::deref(components_[2]);
        }
        throw_runtime_error("Error - VecFieldView not usable");
        return core::deref(components_[0]); // hax
    }

    auto& operator()(Component component) const _PHARE_ALL_FN_ { return getComponent(component); }
    auto& operator()(Component component) _PHARE_ALL_FN_ { return getComponent(component); }


    auto getComponents() const _PHARE_ALL_FN_
    {
        return std::forward_as_tuple((*this)[0], (*this)[1], (*this)[2]);
    }
    auto getComponents() _PHARE_ALL_FN_
    {
        return std::forward_as_tuple((*this)[0], (*this)[1], (*this)[2]);
    }

    auto operator()() const _PHARE_ALL_FN_ { return getComponents(); }
    auto operator()() _PHARE_ALL_FN_ { return getComponents(); }

    auto begin() const { return components_.begin(); }
    auto begin() { return components_.begin(); }

    auto end() const { return components_.end(); }
    auto end() { return components_.end(); }

    auto& operator[](std::size_t i) _PHARE_ALL_FN_ { return core::deref(components_[i]); }
    auto& operator[](std::size_t i) const _PHARE_ALL_FN_ { return core::deref(components_[i]); }

protected:
    Components components_;
};

} // namespace PHARE::core

namespace PHARE
{
namespace core
{
    /** VecField objects represents a 3D vector field in PHARE.
     *
     *  VecField is a ResourcesUser. Like most data objects in PHARE, it does not
     *  own its data, but has pointers to Field resources owned by the SAMRAI system.
     *  Thus VecField has to satisfy the interface required by the ResourcesManager.
     *
     *  VecField class is templated by the type of NdArray Field use and which
     *  physical quantities they represent.
     */

    template<typename NdArrayImpl, typename PhysicalQuantity, typename DataType = double>
    class VecField : public VecFieldView<Field<NdArrayImpl, typename PhysicalQuantity::Scalar>*>
    {
    public:
        using field_type = Field<NdArrayImpl, typename PhysicalQuantity::Scalar>;
        using Super      = VecFieldView<field_type*>;

        using data_view_t = typename NdArrayImpl::view_t;
        using view_t      = VecFieldView<FieldView<data_view_t, typename PhysicalQuantity::Scalar>>;

        using Super::components_;
        using Super::getComponent;
        using Super::getComponents;

        VecField()                                 = delete;
        VecField& Vecfield(VecField const& source) = delete;
        VecField(VecField&& source)                = default;
        VecField& operator=(VecField const& source) = delete;
        VecField& operator=(VecField&& source) = default;

        static constexpr std::size_t dimension = NdArrayImpl::dimension;


        /**
         * @brief builds a VecField from a name and for a specific vector physical quantity
         * @param name is the name of the VecField, components are going to be called name_x, name_y
         * and name_z
         * @param physQty is any of the vector physical quantities available in PHARE
         */
        VecField(std::string const& name, typename PhysicalQuantity::Vector physQty)
            : Super{ConstArray<field_type*, 3>(nullptr)}
            , name_{name}
            , physQties_{PhysicalQuantity::componentsQuantities(physQty)}
        {
        }

        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------

        struct VecFieldProperties
        {
            std::string name;
            typename PhysicalQuantity::Scalar qty;
        };

        using resources_properties = std::vector<VecFieldProperties>;

        resources_properties getFieldNamesAndQuantities() const
        {
            return {{{componentNames_[0], physQties_[0]},
                     {componentNames_[1], physQties_[1]},
                     {componentNames_[2], physQties_[2]}}};
        }

        void setBuffer(std::string const& bufferName, field_type* field)
        {
            if (auto it = nameToIndex_.find(bufferName); it != std::end(nameToIndex_))
                components_[it->second] = field;
            else
                throw std::runtime_error(
                    "VecField Error - invalid component name, cannot set buffer");
        }




        //! return true if the VecField can be used to access component data
        bool isUsable() const
        {
            return std::all_of(std::begin(components_), std::end(components_),
                               [](auto const& c) { return c != nullptr; });
        }

        bool isSettable() const
        {
            return std::all_of(std::begin(components_), std::end(components_),
                               [](auto const& c) { return c == nullptr; });
        }



        void zero()
        {
            if (!isUsable())
                throw std::runtime_error(
                    "Error, cannot zero the VecField because it is not usable");

            for (auto& component : components_)
                component->zero();
        }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        std::string const& name() const { return name_; }


        std::string getComponentName(Component component) const
        {
            switch (component)
            {
                case Component::X: return componentNames_[0];
                case Component::Y: return componentNames_[1];
                case Component::Z: return componentNames_[2];
            }
            throw std::runtime_error("Error - VecField not usable");
        }


        auto components() const
        {
            return std::forward_as_tuple((*this)[0], (*this)[1], (*this)[2]);
        }
        auto components() { return std::forward_as_tuple((*this)[0], (*this)[1], (*this)[2]); }

        auto& operator()(Component component) const { return getComponent(component); }
        auto& operator()(Component component) { return getComponent(component); }

        auto operator()() const { return components(); }
        auto operator()() { return components(); }

        auto& operator[](std::size_t i) { return *components_[i]; }
        auto& operator[](std::size_t i) const { return *components_[i]; }


        void copyData(VecField const& source)
        {
            if (!(isUsable() && source.isUsable()))
                throw std::runtime_error("Error, unusable VecField, cannot copyData");

            for (std::size_t i = 0; i < 3; ++i)
                components_[i]->copyData(*source.components_[i]);
        }


        auto view() { return view_t{*this}; }
        auto view() const { return view_t{*this}; }

    private:
        std::string name_ = "No Name";
        std::array<typename PhysicalQuantity::Scalar, 3> physQties_;

        const std::array<std::string, 3> componentNames_{name_ + "_x", name_ + "_y", name_ + "_z"};

        const std::unordered_map<std::string, std::size_t> nameToIndex_{
            {componentNames_[0], 0u}, {componentNames_[1], 1u}, {componentNames_[2], 2u}};
    };


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
