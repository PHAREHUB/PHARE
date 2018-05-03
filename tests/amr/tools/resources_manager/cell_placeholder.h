#ifndef PHARE_TESTS_AMR_TOOLS_RESSOURCE_CELL_PLACEHOLDER_H
#define PHARE_TESTS_AMR_TOOLS_RESSOURCE_CELL_PLACEHOLDER_H

#include "hybrid/hybrid_quantities.h"
#include "tools/resources_manager.h"

#include <SAMRAI/pdat/CellVariable.h>
#include <string>
#include <utility>


namespace PlaceHolder
{
template<typename T>
class CellVariable : public SAMRAI::pdat::CellVariable<T>
{
public:
    CellVariable(SAMRAI::tbox::Dimension const &dimension, std::string const &name,
                 PHARE::HybridQuantity::Scalar qty, std::string const &layoutType)
        : SAMRAI::pdat::CellVariable<T>(dimension, name)
    {
    }
};

class CellField
{
public:
    CellField(std::string const name, PHARE::HybridQuantity::Scalar hq)
        : name_{name}
        , data_{nullptr}
        , hq_{hq}
    {
    }
    typedef std::vector<std::pair<std::string, PHARE::HybridQuantity::Scalar>>
        resources_properties;

    using field_impl = double;

    resources_properties getFieldNamesAndQuantities() const { return {{name_, hq_}}; }

    void setResources(std::string const resourcesManagerName, field_impl *data)
    {
        if (resourcesManagerName == name_)
        {
            data_ = data;
        }
        else
        {
            std::stringstream stream;
            stream << "Expected name was: " << name_
                   << "Requested name was: " << resourcesManagerName << "\n";
            throw std::runtime_error(stream.str());
        }
    }

    bool isValid() const { return data_ != nullptr; }

private:
    std::string name_;
    field_impl *data_;
    PHARE::HybridQuantity::Scalar hq_;
};
} // namespace PlaceHolder

template<>
struct PHARE::FieldType<PlaceHolder::CellField>
{
    using data_type         = SAMRAI::pdat::CellData<PlaceHolder::CellField::field_impl>;
    using variable_type     = PlaceHolder::CellVariable<PlaceHolder::CellField::field_impl>;
    using internal_type_ptr = PlaceHolder::CellField::field_impl *;
};

template<typename T>
struct ParticlesType
{
};


#endif
