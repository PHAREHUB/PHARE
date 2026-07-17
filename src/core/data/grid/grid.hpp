#ifndef PHARE_CORE_DATA_GRID_GRID_BASE_HPP
#define PHARE_CORE_DATA_GRID_GRID_BASE_HPP


#include "core/def.hpp"
#include "core/data/field/field.hpp"
#include "core/def/phare_config.hpp"

#include <array>
#include <string>
#include <cassert>
#include <optional>
#include <algorithm>


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
public:
    using Super                      = NdArrayImpl;
    auto constexpr static dimension  = NdArrayImpl::dimension;
    auto constexpr static alloc_mode = NdArrayImpl::allocator_mode;

    using value_type             = typename NdArrayImpl::type;
    using physical_quantity_type = PhysicalQuantity;
    using field_type             = Field<dimension, PhysicalQuantity, value_type, alloc_mode>;


    Grid()                              = delete;
    Grid(Grid&& source)                 = default;
    Grid& operator=(Grid&& source)      = delete;
    Grid& operator=(Grid const& source) = delete;

    Grid(std::string const& name, PhysicalQuantity qty,
         std::array<std::uint32_t, dimension> const& dims,
         std::optional<value_type> const value = std::nullopt)
        : Super{dims, value}
        , name_{name}
        , qty_{qty}
    {
    }

    template<typename GridLayout_t>
    Grid(std::string const& name, GridLayout_t const& layout, PhysicalQuantity qty,
         std::optional<value_type> const value = std::nullopt)
        : Grid{name, qty, layout.allocSize(qty), value}
    {
    }

    template<typename... Dims>
    Grid(std::string const& name, PhysicalQuantity qty, Dims... dims)
        requires(std::is_integral_v<Dims> && ...)
        : Grid{name, qty, std::array{dims...}}
    {
        static_assert(sizeof...(Dims) == dimension, "Invalid dimension");
    }


    Grid(Grid const& source) // let field_ default
        : Super{source}
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

    NO_DISCARD operator field_type&() { return field_; }
    NO_DISCARD auto& operator*() { return field_; }
    NO_DISCARD auto& operator*() const { return field_; }

    template<typename, typename>
    friend std::ostream& operator<<(std::ostream& out, Grid const&);

private:
    std::string name_{"No Name"};
    PhysicalQuantity qty_;
    field_type field_{name_, qty_, Super::data(), Super::shape()};
};




template<typename Arr, typename PQ>
struct is_field<Grid<Arr, PQ>> : std::true_type
{
};


template<typename Arr, typename PQ>
inline std::ostream& operator<<(std::ostream& out, Grid<Arr, PQ> const& f)
{
    out << *f;
    return out;
}

template<typename Arr, typename PQ>
inline auto sum_field(Grid<Arr, PQ> const& f)
{
    return sum_field(*f);
}

template<typename Arr, typename PQ>
inline auto sum_not_nan(Grid<Arr, PQ> const& f)
{
    return sum_not_nan(*f);
}


} // namespace PHARE::core

#endif
