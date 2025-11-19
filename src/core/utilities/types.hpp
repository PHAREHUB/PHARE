#ifndef PHARE_CORE_UTILITIES_TYPES_HPP
#define PHARE_CORE_UTILITIES_TYPES_HPP


#include "core/def.hpp"



#include <array>
#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <optional>
#include <algorithm>
#include <stdexcept>




namespace PHARE
{
namespace core
{
    enum class Basis { Magnetic, Cartesian };




    template<typename T>
    NO_DISCARD std::vector<T> arange(T start, T stop, T step = 1)
    {
        std::vector<T> values;
        for (T value = start; value < stop; value += step)
            values.push_back(value);
        return values;
    }

    template<typename T>
    NO_DISCARD T norm(std::array<T, 3> vec)
    {
        auto squarreSum = std::inner_product(std::begin(vec), std::end(vec), std::begin(vec), 0.);
        return std::sqrt(squarreSum);
    }



    enum class Edge { Xmin, Xmax, Ymin, Ymax, Zmin, Zmax };


    template<typename T> // this is so we can use struct {} initialization with
                         // shared_ptrs/forwarding
    struct aggregate_adapter : public T
    {
        template<class... Args>
        aggregate_adapter(Args&&... args)
            : T{std::forward<Args>(args)...}
        {
        }
    };




    template<typename T>
    struct is_std_vector : std::false_type
    {
    };

    template<typename T>
    struct is_std_vector<std::vector<T>> : std::true_type
    {
    };

    template<typename T>
    inline constexpr auto is_std_vector_v = is_std_vector<T>::value;

    template<typename T, std::size_t size>
    struct is_std_array : std::false_type
    {
    };

    template<typename T, std::size_t size>
    struct is_std_array<std::array<T, size>, size> : std::true_type
    {
    };

    template<typename T, std::size_t size>
    inline constexpr auto is_std_array_v = is_std_array<T, size>::value;



    template<typename Tuple, typename Func>
    void apply(Tuple tuple, Func func)
    {
        std::apply([&](auto&... args) { (func(args), ...); }, tuple);
    }

    template<typename Type, std::size_t Size> // std::array::fill is only constexpr in C++20 ffs
    constexpr void fill(Type value, std::array<Type, Size>& array)
    {
        for (std::size_t i = 0; i < Size; i++)
            array[i] = value;
    }

    template<std::size_t Constant>
    class StrongIntegralConstant
    {
    public:
        constexpr decltype(auto) operator()() const { return constant(); }

    protected:
        static constexpr std::integral_constant<std::size_t, Constant> constant{};
    };

    template<std::size_t Constant>
    class DimConst : public StrongIntegralConstant<Constant>
    {
    };
    template<std::size_t Constant>
    class InterpConst : public StrongIntegralConstant<Constant>
    {
    };
    template<std::size_t Constant>
    class RefinedParticlesConst : public StrongIntegralConstant<Constant>
    {
    };

    template<std::size_t To_Size, std::size_t From_Size, typename Type>
    NO_DISCARD constexpr std::array<Type, To_Size>
    sized_array(std::array<Type, From_Size> const& from)
    {
        static_assert(To_Size <= From_Size, "invalid sized_array Size template, too large");

        if constexpr (To_Size == From_Size)
            return from;

        std::array<Type, To_Size> to{};

        for (std::size_t i = 0; i < to.size(); i++)
            to[i] = from[i];

        return to;
    }


    template<std::size_t To_Size, typename... Args>
    NO_DISCARD constexpr auto as_sized_array(Args&&... args)
    {
        auto arr = std::array{std::forward<decltype(args)>(args)...};

        return sized_array<To_Size>(arr);
    }

    template<typename Type, std::size_t size>
    NO_DISCARD constexpr std::array<Type, size> ConstArray(Type val = 0)
    {
        std::array<Type, size> arr{};
        for (uint8_t i = 0; i < size; i++)
            arr[i] = val;
        return arr;
    }

    template<typename Type>
    NO_DISCARD std::vector<Type> displacementFrom(std::vector<Type> const& input)
    {
        std::vector<Type> displs(input.size());
        Type off = 0;
        for (Type i = 0; i < static_cast<Type>(input.size()); i++)
        {
            displs[i] = off;
            off += input[i];
        }
        return displs;
    }

    template<typename T>
    struct StackVar
    {
        using value_type = T;

        T var;
    };

    template<typename T>
    NO_DISCARD std::string to_string_with_precision(T const& a_value, std::size_t const len)
    {
        std::ostringstream out;
        out.precision(len);
        out << std::fixed << a_value;
        auto str = out.str();
        // last digit may be non 0 because of rounding
        // and we know at that decimal it should be so we force it
        str.replace(str.end() - 1, str.end(), 1, '0');
        return out.str();
    }

    template<typename T>
    NO_DISCARD auto to_string_fixed_width(T const& value, std::size_t const& precision,
                                          std::size_t const& width, char const& fill = '0')
    {
        std::ostringstream out;
        out.width(width);
        out.precision(precision);
        out << std::setfill(fill) << std::fixed << value;
        return out.str();
    }

    template<typename T, std::size_t... Is>
    constexpr auto gft_helper(std::index_sequence<Is...> const&&)
        -> decltype(std::make_tuple((Is, std::declval<T>())...));

    template<typename T, std::size_t N>
    constexpr auto get_fixed_tuple() -> decltype(gft_helper<T>(std::make_index_sequence<N>{}));

    template<typename T, std::size_t N>
    using tuple_fixed_type = decltype(get_fixed_tuple<T, N>());



    template<typename T>
    NO_DISCARD T from_string(std::string const& s)
    {
        T t;
        std::stringstream ss(s);
        ss >> t;
        if (ss.fail())
            throw std::runtime_error("Cannot Parse T from: " + s);
        return t;
    }


    NO_DISCARD inline std::optional<std::string> get_env(std::string const& key)
    {
        if (char const* val = std::getenv(key.c_str()))
            return std::string{val};
        return std::nullopt;
    }

    NO_DISCARD inline std::string get_env(std::string const& key, std::string const& _default)
    {
        if (auto e = get_env(key))
            return *e;
        return _default;
    }

    template<typename T>
    NO_DISCARD inline T get_env_as(std::string const& key, T const& t)
    {
        if (auto e = get_env(key))
            return from_string<T>(*e);
        return t;
    }



} // namespace core
} // namespace PHARE


namespace PHARE::core
{
template<typename Container, typename Multiplies = typename Container::value_type>
NO_DISCARD Multiplies product(Container const& container, Multiplies mul = 1)
{
    return std::accumulate(container.begin(), container.end(), mul, std::multiplies<Multiplies>());
}

template<typename Container, typename Return = typename Container::value_type>
NO_DISCARD Return sum(Container const& container, Return r = 0)
{
    return std::accumulate(container.begin(), container.end(), r);
}

template<typename Container, typename F>
NO_DISCARD auto sum_from(Container const& container, F fn)
{
    using value_type  = typename Container::value_type;
    using return_type = std::decay_t<std::invoke_result_t<F const&, value_type const&>>;
    return_type sum   = 0;
    for (auto const& el : container)
        sum += fn(el);
    return sum;
}




template<typename F>
NO_DISCARD auto generate(F&& f, std::size_t from, std::size_t to)
{
    assert(from <= to);
    using value_type = std::decay_t<std::invoke_result_t<F&, std::size_t const&>>;
    std::vector<value_type> v;
    std::size_t count = to - from;
    if (count > 0)
        v.reserve(count);
    for (std::size_t i = from; i < to; ++i)
        v.emplace_back(f(i));
    return v;
}

template<typename F>
NO_DISCARD auto generate(F&& f, std::size_t count)
{
    return generate(std::forward<F>(f), 0, count);
}


template<typename F, typename Container>
NO_DISCARD auto generate(F&& f, Container const& container)
{
    using T          = typename Container::value_type;
    using value_type = std::decay_t<std::invoke_result_t<F&, T&>>;
    std::vector<value_type> v1;
    if (container.size() > 0)
        v1.reserve(container.size());
    for (auto& v : container)
        v1.emplace_back(f(v));
    return v1;
}

template<typename F, typename T>
NO_DISCARD auto generate(F&& f, std::vector<T>&& v)
{
    return generate(std::forward<F>(f), v);
}

template<std::size_t Idx, typename F, typename Type, std::size_t Size>
NO_DISCARD auto constexpr generate_array__(F& f, std::array<Type, Size> const& arr)
{
    return f(arr[Idx]);
}

template<typename Type, std::size_t Size, typename F, std::size_t... Is>
NO_DISCARD auto constexpr generate_array_(F& f, std::array<Type, Size> const& arr,
                                          std::integer_sequence<std::size_t, Is...>)
{
    return std::array{generate_array__<Is>(f, arr)...};
}

template<typename F, typename Type, std::size_t Size>
NO_DISCARD auto constexpr generate(F&& f, std::array<Type, Size> const& arr)
{
    return generate_array_(f, arr, std::make_integer_sequence<std::size_t, Size>{});
}

template<typename T>
auto constexpr all_are(auto&&... ts)
{
    return ((std::is_same_v<T, std::decay_t<decltype(ts)>>) && ...);
}

NO_DISCARD auto constexpr any(auto... bools)
    requires(all_are<bool>(bools...))
{
    return (bools || ...);
}

NO_DISCARD auto constexpr all(auto... bools)
    requires(all_are<bool>(bools...))
{
    return (bools && ...);
}

// calls operator bool() or copies bool
auto constexpr static to_bool = [](auto const& v) { return bool{v}; };


template<typename Container, typename Fn = decltype(to_bool)>
NO_DISCARD auto constexpr all(Container const& container, Fn fn = to_bool)
{
    return std::all_of(container.begin(), container.end(), fn);
}

template<typename Container, typename Fn = decltype(to_bool)>
NO_DISCARD auto constexpr any(Container const& container, Fn fn = to_bool)
{
    return std::any_of(container.begin(), container.end(), fn);
}




template<typename Container, typename Fn = decltype(to_bool)>
NO_DISCARD auto constexpr none(Container const& container, Fn fn = to_bool)
{
    return std::none_of(container.begin(), container.end(), fn);
}


template<typename SignedInt, typename UnsignedInt>
bool diff_sign_int_equals(SignedInt const& i0, UnsignedInt const& i1)
{
    static_assert(std::is_unsigned_v<UnsignedInt>);
    static_assert(std::is_signed_v<SignedInt>);
    static_assert(sizeof(UnsignedInt) >= sizeof(SignedInt), "Bad int comparison!");
    if (i0 < 0)
        return false;
    return static_cast<UnsignedInt>(i0) == i1;
}


template<typename Int0, typename Int1>
bool int_equals(Int0 const& i0, Int1 const& i1)
{
    if constexpr (std::is_same_v<Int0, Int1>)
        return i0 == i1;
    else
    {
        if constexpr (std::is_unsigned_v<Int0> and std::is_signed_v<Int1>)
            return diff_sign_int_equals(i1, i0);
        if constexpr (std::is_unsigned_v<Int1> and std::is_signed_v<Int0>)
            return diff_sign_int_equals(i0, i1);
    }
    // reaching here == compiler error
}



auto inline float_equals(float const& a, float const& b, float diff = 1e-6)
{
    return std::abs(a - b) < diff;
}

auto inline float_equals(double const& a, double const& b, double diff = 1e-12)
{
    return std::abs(a - b) < diff;
}

template<typename T = std::uint16_t>
struct Apply
{
    template<T i>
    auto constexpr operator()()
    {
        return std::integral_constant<T, i>{};
    }
};

template<typename Apply, std::uint16_t... Is>
constexpr auto apply_N(Apply& f, std::integer_sequence<std::uint16_t, Is...> const&)
{
    if constexpr (!std::is_same_v<decltype(f.template operator()<0>()), void>)
        return std::make_tuple(f.template operator()<Is>()...);
    (f.template operator()<Is>(), ...);
}

template<std::uint16_t N, typename Apply>
constexpr auto apply_N(Apply&& f)
{
    return apply_N(f, std::make_integer_sequence<std::uint16_t, N>{});
}

enum class for_N_R_mode {
    make_tuple = 0,
    make_array,
    forward_tuple,
};

template<std::uint16_t N, auto M = for_N_R_mode::make_tuple, typename Fn>
constexpr auto for_N(Fn& fn)
{
    /*  // how to use
        for_N<2>([](auto ic) {
            constexpr auto i = ic();
            // ...
        });
    */

    static_assert(std::is_same_v<decltype(M), for_N_R_mode>);
    using return_type
        = std::decay_t<std::invoke_result_t<Fn, std::integral_constant<std::uint16_t, 0>>>;
    constexpr bool returns = !std::is_same_v<return_type, void>;

    if constexpr (returns)
    {
        return std::apply(
            [&](auto... ics) {
                if constexpr (M == for_N_R_mode::make_tuple)
                    return std::make_tuple(fn(ics)...);
                else if constexpr (M == for_N_R_mode::make_array)
                    return std::array{fn(ics)...};
                else if constexpr (M == for_N_R_mode::forward_tuple)
                    return std::forward_as_tuple(fn(ics)...);
                else
                    throw std::runtime_error("unknown return mode");
            },
            apply_N<N>(Apply{}));
    }
    else
        std::apply([&](auto... ics) { (fn(ics), ...); }, apply_N<N>(Apply{}));
}

template<std::uint16_t N, auto M = for_N_R_mode::make_tuple, typename Fn>
constexpr auto for_N(Fn&& fn)
{
    return for_N<N, M>(fn);
}


template<std::uint16_t N, typename Fn>
constexpr auto for_N_make_array(Fn&& fn)
{
    return for_N<N, for_N_R_mode::make_array>(fn);
}


template<std::uint16_t N, typename Fn>
NO_DISCARD constexpr auto for_N_all(Fn&& fn)
{
    return all(for_N<N, for_N_R_mode::make_array>(fn));
}

template<std::uint16_t N, typename Fn>
NO_DISCARD constexpr auto for_N_any(Fn&& fn)
{
    return any(for_N<N, for_N_R_mode::make_array>(fn));
}



template<typename Tuple, std::size_t... Is>
constexpr auto named_pair_seq_(Tuple, std::index_sequence<Is...> const&&)
    -> decltype(std::make_tuple(
        (Is, std::declval<std::pair<std::string, std::tuple_element_t<Is, Tuple>>>())...));

template<typename... Args>
auto constexpr named_pair_seq()
    -> decltype(named_pair_seq_(std::tuple<Args...>{},
                                std::make_index_sequence<sizeof...(Args)>{}));

template<typename... Args>
using NamedTuple = decltype(named_pair_seq<Args...>());

template<typename... Pairs>
auto make_named_tuple(Pairs&&... pairs)
{
    return std::make_tuple(pairs...);
}



template<typename D>
struct Equals
{
    void operator()(auto& d0) { d = d0; }
    D& d;
};

template<typename D>
struct PlusEquals
{
    void operator()(auto& d0) { d += d0; }
    D& d;
};

} // namespace PHARE::core


#endif // PHARE_CORE_UTILITIES_TYPES_HPP
