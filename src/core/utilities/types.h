#ifndef TYPES_H
#define TYPES_H

#include <array>
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
    using uint32 = std::uint32_t;
    using uint64 = std::uint64_t;
    using int32  = std::int32_t;
    using int64  = std::int64_t;

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

    template<typename T, size_t size>
    struct is_std_array : std::false_type
    {
    };

    template<typename T, size_t size>
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

    template<typename Type, size_t Size> // std::array::fill is only constexpr in C++20 ffs
    constexpr void fill(Type value, std::array<Type, Size>& array)
    {
        for (size_t i = 0; i < Size; i++)
            array[i] = value;
    }

    template<size_t Constant>
    class StrongIntegralConstant
    {
    public:
        constexpr decltype(auto) operator()() const { return constant(); }

    protected:
        static constexpr std::integral_constant<std::size_t, Constant> constant{};
    };

    template<size_t Constant>
    class DimConst : public StrongIntegralConstant<Constant>
    {
    };
    template<size_t Constant>
    class InterpConst : public StrongIntegralConstant<Constant>
    {
    };
    template<size_t Constant>
    class RefinedParticlesConst : public StrongIntegralConstant<Constant>
    {
    };

} // namespace core
} // namespace PHARE

#endif // TYPES_H
