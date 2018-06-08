#ifndef PHARE_FUNCTION_H
#define PHARE_FUNCTION_H

#include <array>
#include <cstddef>
#include <functional>


template<typename R, std::size_t N>
class Function
{
};


template<typename R>
class Function<R, 1>
{
public:
    R operator()(double x) { return f(x); }

private:
    std::function<R(double)> f;
};



template<typename R>
class Function<R, 2>
{
public:
    R operator()(double x, double y) { return f(x, y); }

private:
    std::function<R(double, double)> f;
};


template<typename R>
class Function<R, 3>
{
public:
    R operator()(double x, double y, double z) { return f(x, y, z); }

private:
    std::function<R(double, double, double)> f;
};


template<std::size_t N>
using ScalarFunction = Function<double, N>;

template<std::size_t N>
using VectorFunction = Function<std::array<double, 3>, N>;




#endif
