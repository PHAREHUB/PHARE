#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_H
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_H

#include <array>
#include <utility>

#include "data/field/field.h"
#include "vecfield_component.h"

namespace PHARE
{
//! VecField represents a vector field
/** VecField objects encapsulate 3 Field pointers to the data, one per component.
 *  A VecField object is usable, i.e. can give access to its data, only when
 *  isUsable() returns true. Typical usage of a VecField is to get the components
 *  and work with them.
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

    VecField(std::string const& name, typename PhysicalQuantity::Vector physQty)
        : name_{name}
        , physQties_{PhysicalQuantity::componentsQuantities(physQty)}
        , componentNames_{{name + "_x", name + "_y", name + "_z"}}
    {
    }

    using resources_properties
        = std::vector<std::pair<std::string, typename PhysicalQuantity::Scalar>>;

    using field_type = Field<NdArrayImpl, typename PhysicalQuantity::Scalar>;

    resources_properties getFieldNamesAndQuantities() const
    {
        return {{std::make_pair(componentNames_[0], physQties_[0]),
                 std::make_pair(componentNames_[1], physQties_[1]),
                 std::make_pair(componentNames_[2], physQties_[2])}};
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


    std::string const& name() { return name_; }

    // using field_impl = NdArrayImpl;


    Field<NdArrayImpl, typename PhysicalQuantity::Scalar>& getComponent(Component component)
    {
        switch (component)
        {
            case Component::X: return *xComponent_;
            case Component::Y: return *yComponent_;
            case Component::Z: return *zComponent_;
        }
    }

    Field<NdArrayImpl, typename PhysicalQuantity::Scalar> const&
    getComponent(Component component) const
    {
        switch (component)
        {
            case Component::X: return *xComponent_;
            case Component::Y: return *yComponent_;
            case Component::Z: return *zComponent_;
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

} // namespace PHARE

#endif
