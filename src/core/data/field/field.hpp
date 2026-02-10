#ifndef PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP


#include "core/def.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"


#include <array>
#include <string>
#include <cstddef>
#include <utility>


namespace PHARE::core
{
//! Class Field represents a multidimensional (1,2 or 3D) scalar field
/** Users of Field objects needing to know which physical quantity a specific
 *  Field instance represents can get this info by calling physicalQuantity().
 *  Users may also give a string name to a field object and get a name by calling
 *  name().
 */
template<std::size_t dim, typename PhysicalQuantity, typename Data_t = double>
class Field : public NdArrayView<dim, Data_t>
{
    using Super = NdArrayView<dim, Data_t>;

public:
    auto constexpr static dimension = dim;
    using value_type                = Data_t;
    using physical_quantity_type    = PhysicalQuantity;


    Field(std::string const& name, PhysicalQuantity qty, value_type* data = nullptr,
          std::array<std::uint32_t, dim> const& dims = ConstArray<std::uint32_t, dim>())
        : Super{data, dims}
        , name_{name}
        , qty_{qty}
    {
    }

    Field(Field const& source)            = default;
    Field(Field&& source)                 = default;
    Field& operator=(Field&& source)      = default;
    Field& operator=(Field const& source) = default;

    auto& operator=(Field* src)
    {
        setBuffer(src);
        return *this;
    }


    NO_DISCARD auto& name() const { return name_; }
    NO_DISCARD auto& physicalQuantity() const { return qty_; }

    void copyData(Field const& source) { Super::fill_from(source); }

    void setBuffer(Field* const field)
    {
        auto data = field ? field->data() : nullptr;
        if (data)
        {
            assert(field->name() == this->name());
            Super::setShape(field->shape());
        }
        Super::setBuffer(data);
    }

    bool isUsable() const { return Super::data() != nullptr; }
    bool isSettable() const { return !isUsable(); }


    template<typename... Args>
    NO_DISCARD auto& operator()(Args&&... args)
    {
        PHARE_DEBUG_DO(                                                                 //
            if (!isUsable()) throw std::runtime_error("Field is not usable: " + name_); //
        )
        return super()(std::forward<Args>(args)...);
    }
    template<typename... Args>
    NO_DISCARD auto& operator()(Args&&... args) const
    {
        return const_cast<Field&>(*this)(std::forward<Args>(args)...);
    }

    template<std::size_t, typename, typename>
    friend std::ostream& operator<<(std::ostream& out, Field const&);

private:
    std::string name_{"No Name"};
    PhysicalQuantity qty_;

    Super& super() { return *this; }
    Super const& super() const { return *this; }
};




void print_1d_field(auto& out, auto const& comp)
{
    auto const& shape = comp.shape();

    std::size_t idx = -1;
    for (std::size_t i = 0; i < shape[0]; ++i)
        out << comp.data()[++idx] << ", ";
    out << std::endl;
}

void print_2d_field(auto& out, auto const& comp)
{
    auto const& shape = comp.shape();

    std::size_t idx = -1;
    for (std::size_t i = 0; i < shape[0]; ++i)
    {
        for (std::size_t j = 0; j < shape[1]; ++j)
            out << comp.data()[++idx] << ", ";

        out << std::endl;
    }
    out << std::endl;
}

void print_3d_field(auto& out, auto const& comp)
{
    auto const& shape = comp.shape();

    std::size_t idx = -1;
    for (std::size_t i = 0; i < shape[0]; ++i)
    {
        for (std::size_t j = 0; j < shape[1]; ++j)
        {
            for (std::size_t k = 0; k < shape[2]; ++k)
                out << comp.data()[++idx] << ", ";

            out << std::endl;
        }
        out << std::endl;
    }
    out << std::endl;
}

template<std::size_t dim, typename PQ, typename Data_t>
inline std::ostream& operator<<(std::ostream& out, Field<dim, PQ, Data_t> const& f)
{
    out << f.name() << std::endl;

    if constexpr (dim == 1)
        print_1d_field(out, f);
    if constexpr (dim == 2)
        print_2d_field(out, f);
    if constexpr (dim == 3)
        print_3d_field(out, f);

    return out;
}


} // namespace PHARE::core


#endif
