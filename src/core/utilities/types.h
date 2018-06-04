#ifndef TYPES_H
#define TYPES_H

#include <array>
#include <cinttypes>


namespace PHARE
{
using uint32 = std::uint32_t;
using uint64 = std::uint64_t;
using int32  = std::int32_t;
using int64  = std::int64_t;

enum class Basis { Magnetic, Cartesian };



// using ScalarFunction = double (*)(double x, double y, double z);
// using VectorFunction = std::array<double, 3> (*)(double x, double y, double z);

/**
 * @brief represents a function of 3D spatial coordinates returning a scalar
 */
class ScalarFunction
{
public:
    virtual double operator()(double x, double y, double z) = 0;

    virtual ~ScalarFunction() = default;
};


/**
 * @brief represents a function of 3D spatial coordinates returning a vector
 */
class VectorFunction
{
public:
    virtual std::array<double, 3> operator()(double x, double y, double z) = 0;

    virtual ~VectorFunction() = default;
};



enum class Edge { Xmin, Xmax, Ymin, Ymax, Zmin, Zmax };


} // namespace PHARE

#endif // TYPES_H
