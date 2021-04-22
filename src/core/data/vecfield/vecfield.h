#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_H
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_H

#include <array>
#include <utility>
#include <algorithm>
#include <unordered_map>

#include "core/data/field/field.h"
#include "vecfield_component.h"
#include "core/utilities/meta/meta_utilities.h"

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
    class VecField
    {
    public:
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
            : name_{name}
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

        using field_type = Field<NdArrayImpl, typename PhysicalQuantity::Scalar>;

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
            if (isUsable())
            {
                for (auto& component : components_)
                    component->zero();
            }
            else
            {
                throw std::runtime_error(
                    "Error, cannot zero the VecField because it is not usable");
            }
        }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        std::string const& name() const { return name_; }



        Field<NdArrayImpl, typename PhysicalQuantity::Scalar>& getComponent(Component component)
        {
            if (isUsable())
            {
                switch (component)
                {
                    case Component::X: return *components_[0];
                    case Component::Y: return *components_[1];
                    case Component::Z: return *components_[2];
                }
            }
            throw std::runtime_error("Error - VecField not usable");
        }

        auto getComponents()
        {
            return std::forward_as_tuple(*components_[0], *components_[1], *components_[2]);
        }
        auto getComponents() const
        {
            return std::forward_as_tuple(*components_[0], *components_[1], *components_[2]);
        }


        Field<NdArrayImpl, typename PhysicalQuantity::Scalar> const&
        getComponent(Component component) const
        {
            if (isUsable())
            {
                switch (component)
                {
                    case Component::X: return *components_[0];
                    case Component::Y: return *components_[1];
                    case Component::Z: return *components_[2];
                }
            }
            throw std::runtime_error("Error - VecField not usable");
        }



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



        void copyData(VecField const& source)
        {
            if (isUsable() && source.isUsable())
            {
                for (std::size_t i = 0; i < 3; ++i)
                {
                    components_[i]->copyData(*source.components_[i]);
                }
            }
            else
            {
                throw std::runtime_error("Error, unusable VecField, cannot copyData");
            }
        }


        auto begin() { return std::begin(components_); }

        auto cbegin() const { return std::cbegin(components_); }

        auto end() { return std::end(components_); }

        auto cend() const { return std::cend(components_); }



    private:
        std::string name_ = "No Name";
        std::array<typename PhysicalQuantity::Scalar, 3> physQties_;

        const std::array<std::string, 3> componentNames_{name_ + "_x", name_ + "_y", name_ + "_z"};
        std::array<field_type*, 3> components_{nullptr, nullptr, nullptr};

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
