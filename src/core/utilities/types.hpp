#ifndef TYPES_HPP
#define TYPES_HPP

#include <cassert>
#include <array>
#include <iomanip>
#include <optional>
#include <string>
#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>


#include "cppdict/include/dict.hpp"




namespace PHARE
{
namespace core
{
    enum class Basis { Magnetic, Cartesian };




    template<typename T>
    std::vector<T> arange(T start, T stop, T step = 1)
    {
        std::vector<T> values;
        for (T value = start; value < stop; value += step)
            values.push_back(value);
        return values;
    }

    template<typename T>
    T norm(std::array<T, 3> vec)
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

    template<typename... Args> // this is so we can specialize
    struct type_list           // templates with only the outter most type
    {
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
    constexpr std::array<Type, To_Size> sized_array(std::array<Type, From_Size> const& from)
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
    constexpr auto as_sized_array(Args&&... args)
    {
        auto arr = std::array{std::forward<decltype(args)>(args)...};

        return sized_array<To_Size>(arr);
    }

    template<typename Type, std::size_t size>
    constexpr std::array<Type, size> ConstArray(Type val = 0)
    {
        std::array<Type, size> arr{};
        for (uint8_t i = 0; i < size; i++)
            arr[i] = val;
        return arr;
    }

    template<typename Type>
    std::vector<Type> displacementFrom(std::vector<Type> const& input)
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


        template<typename... Args>
        StackVar(Args&&... args)
            : var(std::forward<Args...>(args...))
        {
        }

        T var;
    };

    template<typename T>
    std::string to_string_with_precision(T const& a_value, std::size_t const len)
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
    auto to_string_fixed_width(T const& value, std::size_t const& precision,
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



    inline std::optional<std::string> get_env(std::string const& key)
    {
        if (const char* val = std::getenv(key.c_str()))
            return std::string{val};
        return std::nullopt;
    }
    inline std::optional<std::string> get_env(std::string&& key) { return get_env(key); }

} // namespace core
} // namespace PHARE


namespace PHARE::core
{
template<typename Container, typename Multiplies = typename Container::value_type>
Multiplies product(Container const& container, Multiplies mul = 1)
{
    return std::accumulate(container.begin(), container.end(), mul, std::multiplies<Multiplies>());
}

template<typename Container, typename F>
auto sum_from(Container const& container, F fn)
{
    using value_type  = typename Container::value_type;
    using return_type = std::decay_t<std::result_of_t<F const&(value_type const&)>>;
    return_type sum   = 0;
    for (auto const& el : container)
        sum += fn(el);
    return sum;
}

template<typename Container, typename Return = typename Container::value_type>
Return sum(Container const& container, Return r = 0)
{
    return std::accumulate(container.begin(), container.end(), r);
}




template<typename F>
auto generate(F&& f, std::size_t from, std::size_t to)
{
    assert(from <= to);
    using value_type = std::decay_t<std::result_of_t<F&(std::size_t const&)>>;
    std::vector<value_type> v;
    std::size_t count = to - from;
    if (count > 0)
        v.reserve(count);
    for (std::size_t i = from; i < to; ++i)
        v.emplace_back(f(i));
    return v;
}

template<typename F>
auto generate(F&& f, std::size_t count)
{
    return generate(std::forward<F>(f), 0, count);
}


template<typename F, typename Container>
auto generate(F&& f, Container& container)
{
    using T          = typename Container::value_type;
    using value_type = std::decay_t<std::result_of_t<F&(T&)>>;
    std::vector<value_type> v1;
    if (container.size() > 0)
        v1.reserve(container.size());
    for (auto& v : container)
        v1.emplace_back(f(v));
    return v1;
}

template<typename F, typename T>
auto generate(F&& f, std::vector<T>&& v)
{
    return generate(std::forward<F>(f), v);
}


template<std::size_t Idx, typename F, typename Type, std::size_t Size>
auto constexpr generate_array__(F& f, std::array<Type, Size>& arr)
{
    return f(arr[Idx]);
}

template<typename Type, std::size_t Size, typename F, std::size_t... Is>
auto constexpr generate_array_(F& f, std::array<Type, Size>& arr,
                               std::integer_sequence<std::size_t, Is...>)
{
    return std::array{generate_array__<Is>(f, arr)...};
}

template<typename F, typename Type, std::size_t Size>
auto constexpr generate(F&& f, std::array<Type, Size> const& arr)
// auto constexpr generate(F&& f, std::array<Type, Size>& arr)
{
    return generate_array_(f, arr, std::make_integer_sequence<std::size_t, Size>{});
}


// calls operator bool() or copies bool
auto constexpr static to_bool = [](auto const& v) { return bool{v}; };


template<typename Container>
auto all(Container const& container)
{
    return std::all_of(container.begin(), container.end(), to_bool);
}

template<typename Container>
auto any(Container const& container)
{
    return std::any_of(container.begin(), container.end(), to_bool);
}

template<typename Container>
auto none(Container const& container)
{
    return std::none_of(container.begin(), container.end(), to_bool);
}



void inline abort_if(bool b)
{
    if (b)
        std::abort();
}


template<typename T = std::size_t>
struct Apply
{
    template<size_t i>
    constexpr auto operator()()
    {
        return std::integral_constant<T, i>{};
    }
};

template<typename Apply, size_t... Is>
constexpr auto apply_N(Apply& f, std::integer_sequence<size_t, Is...> const&)
{
    if constexpr (!std::is_same_v<decltype(f.template operator()<0>()), void>)
        return std::make_tuple(f.template operator()<Is>()...);
    (f.template operator()<Is>(), ...);
}
template<size_t N, typename Apply>
constexpr auto apply_N(Apply&& f)
{
    return apply_N(f, std::make_integer_sequence<size_t, N>{});
}

template<size_t N, typename Fn>
constexpr void for_N(Fn&& fn)
{
    /*
        for_N<2>([](auto ic) {
            constexpr auto i = ic();
            // ...
        });
    */

    std::apply([&](auto... ics) { (fn(ics), ...); }, apply_N<N>(Apply{}));
}




} // namespace PHARE::core


#endif // TYPES_HPP
