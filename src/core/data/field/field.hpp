#ifndef PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP

#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

#include "core/def.hpp"


namespace PHARE::core
{
//! Class Field represents a multidimensional (1,2 or 3D) scalar field
/** Users of Field objects needing to know which physical quantity a specific
 *  Field instance represents can get this info by calling physicalQuantity().
 *  Users may also give a string name to a field object and get a name by calling
 *  name().
 */
template<typename NdArrayImpl, typename PhysicalQuantity>
class Field : public NdArrayImpl
{
    using Super = NdArrayImpl; // only view now
public:
    using value_type             = typename NdArrayImpl::type;
    using physical_quantity_type = PhysicalQuantity;
    using NdArrayImpl::dimension;

    template<std::size_t dim>
    Field(std::string const& name, PhysicalQuantity qty, value_type* data,
          std::array<std::uint32_t, dim> const& dims)
        : Super{data, dims}
        , name_{name}
        , qty_{qty}
    {
    }

    Field(Field const& source)            = default;
    Field(Field&& source)                 = default;
    Field& operator=(Field&& source)      = default;
    Field& operator=(Field const& source) = default;

    NO_DISCARD auto& name() const { return name_; }
    NO_DISCARD auto& physicalQuantity() const { return qty_; }

    void copyData(Field const& source) { Super::fill_from(source); }

private:
    std::string name_{"No Name"};
    PhysicalQuantity qty_;
};




template<typename NdArrayImpl, typename PhysicalQuantity>
void average(Field<NdArrayImpl, PhysicalQuantity> const& f1,
             Field<NdArrayImpl, PhysicalQuantity> const& f2,
             Field<NdArrayImpl, PhysicalQuantity>& avg)
{
    std::transform(std::begin(f1), std::end(f1), std::begin(f2), std::begin(avg),
                   std::plus<double>());

    std::transform(std::begin(avg), std::end(avg), std::begin(avg),
                   [](double x) { return x * 0.5; });
}


} // namespace PHARE::core


#endif
