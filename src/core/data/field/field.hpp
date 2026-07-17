#ifndef PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BASE_HPP

#include "core/def.hpp"
#include "core/def/phare_config.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/ndarray/ndarray_view.hpp"
#include "core/utilities/meta/meta_utilities.hpp"


#include <array>
#include <string>
#include <cstddef>
// #include <stdexcept>


namespace PHARE::core
{
template<typename T, typename Attempt = void>
struct has_physicalQuantity : std::false_type
{
};

template<typename T>
struct has_physicalQuantity<T,
                            core::tryToInstanciate<decltype(std::declval<T>().physicalQuantity())>>
    : std::true_type
{
};
template<typename T>
auto constexpr has_physicalQuantity_v = has_physicalQuantity<T>::value;


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

    NO_DISCARD auto& physicalQuantity() const _PHARE_ALL_FN_ { return qty_; }

    bool isUsable() const { return Super::data() != nullptr; }
    bool isSettable() const { return !isUsable(); }

    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }

    Super& super() _PHARE_ALL_FN_ { return *this; }
    Super const& super() const _PHARE_ALL_FN_ { return *this; }

    template<auto>
    friend std::ostream& operator<<(std::ostream& out, Field const&);


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

public:
    using Super = basic::Field<FieldOpts<PhysicalQuantity, Data_t>{dim, alloc_mode_}>;
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


    NO_DISCARD auto& name() const { return name_; }

    void copyData(Field const& source) { Super::fill_from(source); }

    void setBuffer(std::nullptr_t ptr) _PHARE_ALL_FN_ { setBuffer(static_cast<Field*>(nullptr)); }

    template<typename FieldLike>
    void setBuffer(FieldLike* const field) _PHARE_ALL_FN_
    {
        auto data = field ? field->data() : nullptr;
        if (data)
        {
            assert(field->name() == this->name());
            Super::setShape(field->shape());
        }
        Super::setBuffer(data);
    }

    void setData(Data_t* const data) _PHARE_ALL_FN_ { Super::setBuffer(data); }

    bool isUsable() const { return Super::data() != nullptr; }
    bool isSettable() const { return !isUsable(); }


    template<typename... Args>
    NO_DISCARD auto& operator()(Args&&... args) _PHARE_ALL_FN_
    {
        if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            PHARE_DEBUG_DO(                                                                 //
                if (!isUsable()) throw std::runtime_error("Field is not usable: " + name_); //
            )
        }

        return super()(args...);
    }
    template<typename... Args>
    NO_DISCARD auto const& operator()(Args&&... args) const _PHARE_ALL_FN_
    {
        return super()(args...);
    }

    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }



private:
    std::string name_{"No Name"};

    Super& super() _PHARE_ALL_FN_ { return *this; }
    Super const& super() const _PHARE_ALL_FN_ { return *this; }
};


template<typename FieldLike_t, typename Data_t>
auto make_field_from(FieldLike_t const& field, Data_t* data)
{
    using Field_t
        = Field<FieldLike_t::dimension, typename FieldLike_t::physical_quantity_type, Data_t>;

    return Field_t{field.name(), field.physicalQuantity(), data, field.shape()};
}


template<typename T>
struct is_field : std::false_type
{
};

template<std::size_t dim, typename D>
struct is_field<NdArrayView<dim, D>> : std::true_type
{
};

template<std::size_t dim, typename PQ, typename D, auto am>
struct is_field<Field<dim, PQ, D, am>> : std::true_type
{
};
template<auto opts>
struct is_field<basic::Field<opts>> : std::true_type
{
};


template<typename T>
auto static constexpr is_field_v = is_field<T>::value;




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



template<auto opts>
inline std::ostream& operator<<(std::ostream& out, basic::Field<opts> const& f)
{
    // out << f.name() << std::endl;

    if constexpr (opts.dimension == 1)
        print_1d_field(out, f);
    if constexpr (opts.dimension == 2)
        print_2d_field(out, f);
    if constexpr (opts.dimension == 3)
        print_3d_field(out, f);

    return out;
}

template<auto opts>
inline auto sum_field(basic::Field<opts> const& f)
{
    return sum(f);
}

template<auto opts>
inline auto sum_not_nan(basic::Field<opts> const& f)
{
    typename basic::Field<opts>::value_type s = 0;
    for (auto const& v : f)
        if (!std::isnan(v))
            s += v;
    return s;
}


template<auto opts>
void check_field(basic::Field<opts> const& f, auto const& layout)
{
#if PHARE_DEBUG
    for (auto const& bix : layout.domainBoxFor(f))
    {
        // if (std::isnan(f(bix)))
        // {
        //     PHARE_LOG_LINE_SS("\n" << f);
        // }
        // assert(not std::isnan(f(bix)));
    }
#endif // PHARE_DEBUG
}


template<auto opts>
void check_field(basic::Field<opts> const& f)
{
#if PHARE_DEBUG
    for (auto const& v : f)
    {
        if (std::isnan(v))
        {
            PHARE_LOG_LINE_SS("\n" << f);
        }
        assert(not std::isnan(v));
        // assert(v == 0 or std::abs(v) > 1e-40);
    }
#endif // PHARE_DEBUG
}


template<typename Field_t>
void check_field(std::vector<Field_t> const& vec)
{
#if PHARE_DEBUG
    for (auto const& f : vec)
        check_field(deref(f));
#endif // PHARE_DEBUG
}


} // namespace PHARE::core


#endif
