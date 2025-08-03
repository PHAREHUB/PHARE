#ifndef PHARE_CORE_DATA_GRID_GRID_BASE_HPP
#define PHARE_CORE_DATA_GRID_GRID_BASE_HPP


#include "core/def.hpp"
#include "core/data/field/field.hpp"


#include <array>
#include <string>
#include <cstddef>
#include <cassert>

namespace PHARE::core
{


/* Grid is the structure owning the field type memory via its inheritance from NdArrayImpl
Grid exists to decouple the usage of memory by computing routines from the allocation of
memory. Components needing to own/allocate memory will use a Grid.
On the contrary, components that just need to manipulate data (and may not be able to manipulate
objects encapsulating allocating objects such as vectors) will access it through a Field view. For
convenience, Grid can spawn its own Field view.
*/
template<typename NdArrayImpl, typename PhysicalQuantity>
class Grid : public NdArrayImpl
{
    using Super = NdArrayImpl;

public:
    using value_type             = typename NdArrayImpl::type;
    using physical_quantity_type = PhysicalQuantity;
    using NdArrayImpl::dimension;
    using field_type = Field<dimension, PhysicalQuantity, value_type>;


    Grid()                              = delete;
    Grid(Grid&& source)                 = default;
    Grid& operator=(Grid&& source)      = delete;
    Grid& operator=(Grid const& source) = delete;

    template<typename... Dims>
    Grid(std::string const& name, PhysicalQuantity qty, Dims... dims)
        : Super{dims...}
        , name_{name}
        , qty_{qty}
    {
        static_assert(sizeof...(Dims) == dimension, "Invalid dimension");
    }

    template<FloatingPoint U = value_type, std::size_t dim>
    Grid(std::string const& name, PhysicalQuantity qty, std::array<std::uint32_t, dim> const& dims,
         value_type value = static_cast<U>(std::nan("")))
        : Super{dims, value}
        , name_{name}
        , qty_{qty}
    {
    }

    template<FloatingPoint U = value_type, typename GridLayout_t>
    Grid(std::string const& name, GridLayout_t const& layout, PhysicalQuantity qty,
         value_type value = static_cast<U>(std::nan("")))
        : Super{layout.allocSize(qty), value}
        , name_{name}
        , qty_{qty}
    {
    }

    template<std::size_t dim>
        requires(!FloatingPoint<value_type>)
    Grid(std::string const& name, PhysicalQuantity qty, std::array<std::uint32_t, dim> const& dims)
        : Super{dims}
        , name_{name}
        , qty_{qty}
    {
    }

    template<typename GridLayout_t>
        requires(!FloatingPoint<value_type>)
    Grid(std::string const& name, GridLayout_t const& layout, PhysicalQuantity qty)
        : Super{layout.allocSize(qty)}
        , name_{name}
        , qty_{qty}
    {
    }
    Grid(Grid const& source) // let field_ default
        : Super{source.shape()}
        , name_{source.name()}
        , qty_{source.physicalQuantity()}
    {
    }

    NO_DISCARD std::string name() const { return name_; }

    NO_DISCARD constexpr PhysicalQuantity physicalQuantity() const { return qty_; }

    template<typename That>
    void copyData(That const& that)
    {
        if (for_N_any<dimension>([&](auto i) { return this->shape()[i] != that.shape()[i]; }))
            throw std::runtime_error("Grid::copyData: Incompatible input shape");
        std::copy(that.data(), that.data() + Super::size(), Super::data());
    }

    void zero() { field_.zero(); } // is always usable

    // returns view when getting address of this object, could be misleading, but convenient
    NO_DISCARD auto operator&() { return &field_; }
    NO_DISCARD auto operator&() const { return &field_; }

private:
    std::string name_{"No Name"};
    PhysicalQuantity qty_;
    field_type field_{name_, qty_, Super::data(), Super::shape()};
};



} // namespace PHARE::core

#endif
