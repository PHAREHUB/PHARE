#ifndef PHARE_CORE_DATA_GRID_GRID_BASE_HPP
#define PHARE_CORE_DATA_GRID_GRID_BASE_HPP

#include <array>
#include <cstddef>
#include <cassert>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

#include "core/def.hpp"
#include "core/data/field/field.hpp"

namespace PHARE::core
{


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
    Grid(Grid const& source)            = delete;
    Grid(Grid&& source)                 = default;
    Grid& operator=(Grid&& source)      = delete;
    Grid& operator=(Grid const& source) = delete;

    template<typename... Dims>
    Grid(std::string const& name, PhysicalQuantity qty, Dims... dims)
        : Super{dims...}
        , name_{name}
        , qty_{qty}
        , field_{name, qty, Super::data(), Super::shape()}
    {
        static_assert(sizeof...(Dims) == dimension, "Invalid dimension");
    }

    template<std::size_t dim>
    Grid(std::string const& name, PhysicalQuantity qty, std::array<std::uint32_t, dim> const& dims)
        : Super{dims}
        , name_{name}
        , qty_{qty}
        , field_{name, qty, Super::data(), Super::shape()}
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

    // returns view when getting address of this object, could be misleading, but convenient
    NO_DISCARD auto operator&() { return &field_; }
    NO_DISCARD auto operator&() const { return &field_; }

private:
    std::string name_{"No Name"};
    PhysicalQuantity qty_;
    field_type field_;
};




template<typename NdArrayImpl, typename PhysicalQuantity>
void average(Grid<NdArrayImpl, PhysicalQuantity> const& f1,
             Grid<NdArrayImpl, PhysicalQuantity> const& f2,
             Grid<NdArrayImpl, PhysicalQuantity>& avg)
{
    std::transform(std::begin(f1), std::end(f1), std::begin(f2), std::begin(avg),
                   std::plus<double>());

    std::transform(std::begin(avg), std::end(avg), std::begin(avg),
                   [](double x) { return x * 0.5; });
}



} // namespace PHARE::core

#endif
