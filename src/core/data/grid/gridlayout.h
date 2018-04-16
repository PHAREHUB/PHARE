#ifndef PHARE_CORE_DATA_GRID_GRIDLAYOUT_H
#define PHARE_CORE_DATA_GRID_GRIDLAYOUT_H

#include <array>
#include <cstddef>

#include <hybrid/hybrid_quantities.h>

namespace PHARE
{
enum class Centering { primal, dual };

enum LayoutType { Yee, Normal };


constexpr int centering2int(Centering c)
{
    return static_cast<int>(c);
}

// mock for a GridLayout
// here all that matters is that GridLayout is able to give me
// at compile time the centering of a given HybridQuantity
template<LayoutType, std::size_t dim>
class GridLayout
{
public:
    static constexpr std::array<Centering, dim> centering(HybridQuantity::Quantity qty);
};


// this function returns the centering of a given hybrid quantity
// in all three directions. A priori this should be inlined/optimized
// at compile time
template<>
constexpr std::array<Centering, 3> GridLayout<Yee, 3>::centering(HybridQuantity::Quantity qty)
{
    switch (qty)
    {
        case HybridQuantity::Quantity::Ex:
            return {{Centering::dual, Centering::primal, Centering::primal}};
        case HybridQuantity::Quantity::Ey:
            return {{Centering::primal, Centering::dual, Centering::primal}};
        case HybridQuantity::Quantity::Ez:
            return {{Centering::primal, Centering::primal, Centering::dual}};
        case HybridQuantity::Quantity::Bx:
            return {{Centering::primal, Centering::dual, Centering::dual}};
        case HybridQuantity::Quantity::By:
            return {{Centering::dual, Centering::primal, Centering::dual}};
        case HybridQuantity::Quantity::Bz:
            return {{Centering::dual, Centering::dual, Centering::primal}};
    }
}


template<>
constexpr std::array<Centering, 2> GridLayout<Yee, 2>::centering(HybridQuantity::Quantity qty)
{
    switch (qty)
    {
        case HybridQuantity::Quantity::Ex: return {Centering::dual, Centering::primal};
        case HybridQuantity::Quantity::Ey: return {Centering::primal, Centering::dual};
        case HybridQuantity::Quantity::Ez: return {Centering::primal, Centering::primal};
        case HybridQuantity::Quantity::Bx: return {Centering::primal, Centering::dual};
        case HybridQuantity::Quantity::By: return {Centering::dual, Centering::primal};
        case HybridQuantity::Quantity::Bz: return {Centering::dual, Centering::dual};
    }
}


template<>
constexpr std::array<Centering, 1> GridLayout<Yee, 1>::centering(HybridQuantity::Quantity qty)
{
    switch (qty)
    {
        case HybridQuantity::Quantity::Ex: return {Centering::dual};
        case HybridQuantity::Quantity::Ey: return {Centering::primal};
        case HybridQuantity::Quantity::Ez: return {Centering::primal};
        case HybridQuantity::Quantity::Bx: return {Centering::primal};
        case HybridQuantity::Quantity::By: return {Centering::dual};
        case HybridQuantity::Quantity::Bz: return {Centering::dual};
    }
}


} // namespace PHARE

#endif
