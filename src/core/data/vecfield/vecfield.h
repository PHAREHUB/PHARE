#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_H
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_H

#include <array>
#include <utility>

#include "data/field/field.h"
#include "vecfield_component.h"

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



        /**
         * @brief builds a VecField from a name and for a specific vector physical quantity
         * @param name is the name of the VecField, components are going to be called name_x, name_y
         * and name_z
         * @param physQty is any of the vector physical quantities available in PHARE
         */
        VecField(std::string const& name, typename PhysicalQuantity::Vector physQty)
            : name_{name}
            , physQties_{PhysicalQuantity::componentsQuantities(physQty)}
            , componentNames_{{name + "_x", name + "_y", name + "_z"}}
        {
        }

        static constexpr std::size_t dimension = NdArrayImpl::dimension;

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
            if (bufferName == componentNames_[0])
            {
                xComponent_ = field;
            }
            else if (bufferName == componentNames_[1])
            {
                yComponent_ = field;
            }
            else if (bufferName == componentNames_[2])
            {
                zComponent_ = field;
            }
        }

        //! return true if the VecField can be used to access component data
        bool isUsable() const
        {
            return xComponent_ != nullptr && yComponent_ != nullptr && zComponent_ != nullptr;
        }

        bool isSettable() const
        {
            return xComponent_ == nullptr && yComponent_ == nullptr && zComponent_ == nullptr;
        }



        void zero()
        {
            if (isUsable())
            {
                xComponent_->zero();
                yComponent_->zero();
                zComponent_->zero();
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
                    case Component::X: return *xComponent_;
                    case Component::Y: return *yComponent_;
                    case Component::Z: return *zComponent_;
                }
            }
            throw std::runtime_error("Error - VecField not usable");
        }




        Field<NdArrayImpl, typename PhysicalQuantity::Scalar> const&
        getComponent(Component component) const
        {
            if (isUsable())
            {
                switch (component)
                {
                    case Component::X: return *xComponent_;
                    case Component::Y: return *yComponent_;
                    case Component::Z: return *zComponent_;
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
                xComponent_->copyData(*source.xComponent_);
                yComponent_->copyData(*source.yComponent_);
                zComponent_->copyData(*source.zComponent_);
            }
            else
            {
                throw std::runtime_error("Error, unusable VecField, cannot copyData");
            }
        }


    private:
        std::string name_ = "No Name";
        std::array<typename PhysicalQuantity::Scalar, 3> physQties_;
        std::array<std::string, 3> componentNames_;
        field_type* xComponent_ = nullptr;
        field_type* yComponent_ = nullptr;
        field_type* zComponent_ = nullptr;
    };
} // namespace core
} // namespace PHARE

#endif
