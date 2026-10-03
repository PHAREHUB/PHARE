#ifndef PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP

#include "core/def.hpp"
#include "core/data/vector.hpp"
#include "core/data/ndarray/ndarray_view.hpp"

#include <array>
#include <string>
#include <cstddef>

namespace PHARE::core
{
template<typename T>
concept HasPhysicalQuantity = requires(T t) { t.physicalQuantity(); };

template<typename T>
auto constexpr has_physicalQuantity_v = HasPhysicalQuantity<T>;

template<typename PhysicalQuantity, typename Data_t = double>
struct FieldOpts
{
    using physical_quantity_type = PhysicalQuantity;
    using value_type             = Data_t;

    std::size_t dimension;
    AllocatorMode alloc_mode = AllocatorMode::CPU;
};

} // namespace PHARE::core

namespace PHARE::core::basic
{

template<auto opts> // NO STRINGS OR STD LIB!
class Field : public NdArrayView<opts.dimension, typename decltype(opts)::value_type>
{
    using FieldOpts = decltype(opts);

public:
    using physical_quantity_type     = FieldOpts::physical_quantity_type;
    using value_type                 = typename FieldOpts::value_type;
    using Super                      = NdArrayView<opts.dimension, value_type>;
    auto constexpr static alloc_mode = opts.alloc_mode;
    auto constexpr static dimension  = opts.dimension;

    Field(physical_quantity_type qty, value_type* data = nullptr,
          std::array<std::uint32_t, opts.dimension> const& dims
          = ConstArray<std::uint32_t, opts.dimension>())
        : Super{data, dims}
        , qty_{qty}
    {
    }

    NO_DISCARD auto& physicalQuantity() const { return qty_; }

    bool isUsable() const { return Super::data() != nullptr; }
    bool isSettable() const { return !isUsable(); }

    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }

    Super& super() { return *this; }
    Super const& super() const { return *this; }

protected:
    physical_quantity_type qty_;
};

} // namespace PHARE::core::basic

namespace PHARE::core
{
//! Class Field represents a multidimensional (1,2 or 3D) scalar field
/** Users of Field objects needing to know which physical quantity a specific
 *  Field instance represents can get this info by calling physicalQuantity().
 *  Users may also give a string name to a field object and get a name by calling
 *  name().
 */
template<std::size_t dim, typename PhysicalQuantity, typename Data_t = double,
         auto alloc_mode_ = AllocatorMode::CPU>
class Field : public basic::Field<FieldOpts<PhysicalQuantity, Data_t>{dim, alloc_mode_}>
{
    static_assert(std::is_same_v<decltype(alloc_mode_), AllocatorMode>);
    using Super = basic::Field<FieldOpts<PhysicalQuantity, Data_t>{dim, alloc_mode_}>;

public:
    auto constexpr static dimension  = dim;
    auto constexpr static alloc_mode = alloc_mode_;
    using value_type                 = Data_t;
    using physical_quantity_type     = PhysicalQuantity;

    Field(std::string const& name, PhysicalQuantity qty, value_type* data = nullptr,
          std::array<std::uint32_t, dim> const& dims = ConstArray<std::uint32_t, dim>())
        : Super{qty, data, dims}
        , name_{name}
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

    template<typename FieldLike>
    void setBuffer(FieldLike* const field)
    {
        auto data = field ? field->data() : nullptr;
        if (data)
        {
            assert(field->name() == this->name());
            Super::setShape(field->shape());
        }
        Super::setBuffer(data);
    }

    void copyData(Field const& source) { Super::fill_from(source); }

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

    void setBuffer(std::nullptr_t ptr) { setBuffer(static_cast<Field*>(nullptr)); }
    void setData(Data_t* const data) { Super::setBuffer(data); }

    NO_DISCARD auto& name() const { return name_; }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }

private:
    std::string name_{"No Name"};

    Super& super() { return *this; }
    Super const& super() const { return *this; }
};

template<typename FieldLike_t, typename Data_t>
auto make_field_from(FieldLike_t const& field, Data_t* data)
{
    using Field_t
        = Field<FieldLike_t::dimension, typename FieldLike_t::physical_quantity_type, Data_t>;

    return Field_t{field.name(), field.physicalQuantity(), data, field.shape()};
}

template<typename T>
inline constexpr bool is_field_v = is_ndarray_v<T>;

} // namespace PHARE::core

#endif
