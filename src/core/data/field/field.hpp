
#ifndef PHARE_CORE_DATA_FIELD_FIELD_BASE_H
#define PHARE_CORE_DATA_FIELD_FIELD_BASE_H

#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

#include "core/data/ndarray/ndarray_vector.hpp"




namespace PHARE
{
namespace core
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
    public:
        using impl_type              = NdArrayImpl;
        using type                   = typename NdArrayImpl::type;
        using physical_quantity_type = PhysicalQuantity;
        using NdArrayImpl::dimension;


        Field()                    = delete;
        Field(Field const& source) = delete;
        Field(Field&& source)      = default;
        Field& operator=(Field&& source) = delete;
        Field& operator=(Field const& source) = delete;

        template<typename... Dims>
        Field(std::string name, PhysicalQuantity qty, Dims... dims)
            : NdArrayImpl{dims...}
            , name_{std::move(name)}
            , qty_{qty}
        {
            static_assert(sizeof...(Dims) == NdArrayImpl::dimension, "Invalid dimension");
        }

        template<std::size_t dim>
        Field(std::string name, PhysicalQuantity qty, std::array<std::uint32_t, dim> const& dims)
            : NdArrayImpl{dims}
            , name_{std::move(name)}
            , qty_{qty}
        {
        }

        std::string name() const { return name_; }


        constexpr PhysicalQuantity physicalQuantity() const { return qty_; }

        void copyData(Field const& source)
        {
            static_cast<NdArrayImpl&>(*this) = static_cast<NdArrayImpl const&>(source);
        }


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


} // namespace core
} // namespace PHARE
#endif
