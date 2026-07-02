#ifndef PHARE_CORE_DATA_PATCH_FIELD_ACCESSOR_HPP
#define PHARE_CORE_DATA_PATCH_FIELD_ACCESSOR_HPP

#include "core/data/vecfield/vecfield.hpp"

namespace PHARE::core
{

/**
 * @brief Abstract interface for accessing fields on a patch by physical quantity.
 *
 * Provides a model-agnostic way to retrieve scalar and vector fields from a patch
 * at runtime, given their physical quantity enum value. Concrete implementations
 * live in the amr layer where patch data storage (e.g. SAMRAI) is available;
 * the core layer only sees this interface.
 *
 * @tparam FieldT The scalar field type.
 * @tparam PhysicalQuantityT The physical quantity category (e.g. MHDQuantity).
 */
template<typename FieldT, typename PhysicalQuantityT>
class IPatchFieldAccessor
{
public:
    using scalar_quantity_type = PhysicalQuantityT::Scalar;
    using vector_quantity_type = PhysicalQuantityT::Vector;
    using field_type           = FieldT;
    using vectorfield_type     = VecField<FieldT, PhysicalQuantityT>;

    virtual ~IPatchFieldAccessor()                                       = default;
    virtual field_type& getField(scalar_quantity_type qty) const         = 0;
    virtual vectorfield_type getVecField(vector_quantity_type qty) const = 0;
};

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_PATCH_FIELD_ACCESSOR_HPP
