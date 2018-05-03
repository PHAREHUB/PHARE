#ifndef PHARE_CORE_DATA_GRID_GRIDLAYOUT_H
#define PHARE_CORE_DATA_GRID_GRIDLAYOUT_H

#include <array>
#include <cstddef>

#include <hybrid/hybrid_quantities.h>

namespace PHARE
{
enum class QtyCentering { primal, dual };

enum Layout { Yee, Normal };


constexpr int centering2int(QtyCentering c)
{
    return static_cast<int>(c);
}

// mock for a GridLayout
// here all that matters is that GridLayout is able to give me
// at compile time the centering of a given HybridQuantity
template<Layout, std::size_t dim>
class GridLayout
{
public:
    static constexpr std::array<QtyCentering, dim> centering(HybridQuantity::Scalar qty);
};


// this function returns the centering of a given hybrid quantity
// in all three directions. A priori this should be inlined/optimized
// at compile time
template<>
constexpr std::array<QtyCentering, 3> GridLayout<Yee, 3>::centering(HybridQuantity::Scalar qty)
{
    switch (qty)
    {
        case HybridQuantity::Scalar::Ex:
            return {{QtyCentering::dual, QtyCentering::primal, QtyCentering::primal}};
        case HybridQuantity::Scalar::Ey:
            return {{QtyCentering::primal, QtyCentering::dual, QtyCentering::primal}};
        case HybridQuantity::Scalar::Ez:
            return {{QtyCentering::primal, QtyCentering::primal, QtyCentering::dual}};
        case HybridQuantity::Scalar::Bx:
            return {{QtyCentering::primal, QtyCentering::dual, QtyCentering::dual}};
        case HybridQuantity::Scalar::By:
            return {{QtyCentering::dual, QtyCentering::primal, QtyCentering::dual}};
        case HybridQuantity::Scalar::Bz:
            return {{QtyCentering::dual, QtyCentering::dual, QtyCentering::primal}};
    }
}


template<>
constexpr std::array<QtyCentering, 2> GridLayout<Yee, 2>::centering(HybridQuantity::Scalar qty)
{
    switch (qty)
    {
        case HybridQuantity::Scalar::Ex: return {QtyCentering::dual, QtyCentering::primal};
        case HybridQuantity::Scalar::Ey: return {QtyCentering::primal, QtyCentering::dual};
        case HybridQuantity::Scalar::Ez: return {QtyCentering::primal, QtyCentering::primal};
        case HybridQuantity::Scalar::Bx: return {QtyCentering::primal, QtyCentering::dual};
        case HybridQuantity::Scalar::By: return {QtyCentering::dual, QtyCentering::primal};
        case HybridQuantity::Scalar::Bz: return {QtyCentering::dual, QtyCentering::dual};
    }
}


template<>
constexpr std::array<QtyCentering, 1> GridLayout<Yee, 1>::centering(HybridQuantity::Scalar qty)
{
    switch (qty)
    {
        case HybridQuantity::Scalar::Ex: return {QtyCentering::dual};
        case HybridQuantity::Scalar::Ey: return {QtyCentering::primal};
        case HybridQuantity::Scalar::Ez: return {QtyCentering::primal};
        case HybridQuantity::Scalar::Bx: return {QtyCentering::primal};
        case HybridQuantity::Scalar::By: return {QtyCentering::dual};
        case HybridQuantity::Scalar::Bz: return {QtyCentering::dual};
    }
}


} // namespace PHARE

#endif
