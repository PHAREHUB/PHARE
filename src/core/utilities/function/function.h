#ifndef PHARE_FUNCTION_H
#define PHARE_FUNCTION_H

#include <array>
#include <cstddef>
#include <functional>


namespace PHARE
{
namespace core
{
    template<typename R, std::size_t N>
    class Function
    {
    };


    template<typename R>
    class Function<R, 1>
    {
    public:
        explicit Function(std::function<R(double)> func)
            : f{func}
        {
        }
        Function() = default;

        Function(Function<R, 1> const& other) = default;

        R operator()(double x) const { return f(x); }

    private:
        std::function<R(double)> f;
    };



    template<typename R>
    class Function<R, 2>
    {
    public:
        explicit Function(std::function<R(double, double)> func)
            : f{func}
        {
        }

        Function() = default;

        Function(Function<R, 2> const& other) = default;

        R operator()(double x, double y) const { return f(x, y); }

    private:
        std::function<R(double, double)> f;
    };


    template<typename R>
    class Function<R, 3>
    {
    public:
        explicit Function(std::function<R(double, double, double)> func)
            : f{func}
        {
        }

        Function() = default;

        Function(Function<R, 3> const& other) = default;

        R operator()(double x, double y, double z) const { return f(x, y, z); }

    private:
        std::function<R(double, double, double)> f;
    };


    template<std::size_t N>
    using ScalarFunction = Function<double, N>;

    template<std::size_t N>
    using VectorFunction = Function<std::array<double, 3>, N>;

} // namespace core
} // namespace PHARE



#endif
