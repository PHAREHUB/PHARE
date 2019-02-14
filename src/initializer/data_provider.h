#ifndef DATA_PROVIDER_H
#define DATA_PROVIDER_H

#include "cppdict/include/dict.hpp"
//#include "models/physical_state.h"

#include <map>
#include <string>
#include <variant>

namespace PHARE
{
namespace initializer
{
    template<typename ReturnType, std::size_t dim>
    struct ScalarFunctionHelper
    {
    };

    template<>
    struct ScalarFunctionHelper<double, 1>
    {
        using type = std::function<double(double)>;
    };

    template<>
    struct ScalarFunctionHelper<double, 2>
    {
        using type = std::function<double(double)>;
    };

    template<>
    struct ScalarFunctionHelper<double, 3>
    {
        using type = std::function<double(double, double, double)>;
    };

    template<typename ReturnType, std::size_t dim>
    struct VectorFunctionHelper
    {
    };

    template<>
    struct VectorFunctionHelper<std::array<double, 3>, 1>
    {
        using type = std::function<std::array<double, 3>(double)>;
    };

    template<>
    struct VectorFunctionHelper<std::array<double, 3>, 2>
    {
        using type = std::function<std::array<double, 3>(double, double)>;
    };

    template<>
    struct VectorFunctionHelper<std::array<double, 3>, 3>
    {
        using type = std::function<std::array<double, 3>(double, double, double)>;
    };

    template<std::size_t dim>
    using ScalarFunction = typename ScalarFunctionHelper<double, dim>::type;

    template<std::size_t dim>
    using VectorFunction = typename VectorFunctionHelper<std::array<double, 3>, dim>::type;



    template<std::size_t dim>
    using PHAREDict = cppdict::Dict<int, double, std::size_t, std::string, ScalarFunction<dim>,
                                    VectorFunction<dim>>;



    extern PHAREDict<1> phareDict1D;
    extern PHAREDict<2> phareDict2D;
    extern PHAREDict<3> phareDict3D;

    template<std::size_t dim>
    auto& dict()
    {
        static_assert(dim >= 1 and dim <= 3, "error, invalid dimension in dict<dim>()");
        if constexpr (dim == 1)
        {
            return phareDict1D;
        }
        else if constexpr (dim == 2)
        {
            return phareDict2D;
        }
        else if constexpr (dim == 3)
        {
            return phareDict3D;
        }
    }

    class DataProvider
    {
    public:
        /**
         * @brief read will read the data from whatever source specific DataProvider classes
         * implement. readData() must be called before passing
         */
        virtual void read()     = 0;
        virtual ~DataProvider() = default;
    };

} // namespace initializer

} // namespace PHARE

#endif // DATA_PROVIDER_H
