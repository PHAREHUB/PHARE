#ifndef PHARE_CORE_GRID_GRIDLAYOUTDEFS_H
#define PHARE_CORE_GRID_GRIDLAYOUTDEFS_H


#include "hybrid/hybrid_quantities.h"
#include "utilities/types.h"

namespace PHARE
{
enum class Direction { X, Y, Z };

enum class QtyCentering { primal, dual };


struct WeightPoint
{
    int32 ix, iy, iz;
    double coef;
};


using LinearCombination = std::vector<WeightPoint>;

enum class Layout { Yee };

/**
 * @brief gridDataT provides constants used to initialize:
 * - hybridQuantity centerings
 * - physical start/end indexes
 * - ghost start/end indexes
 * - numbers of padding cells and physical cells
 */
struct gridDataT
{
    static constexpr Direction dirX = Direction::X;
    static constexpr Direction dirY = Direction::Y;
    static constexpr Direction dirZ = Direction::Z;

    static constexpr QtyCentering primal = QtyCentering::primal;
    static constexpr QtyCentering dual   = QtyCentering::dual;

    static constexpr uint32 idirX = static_cast<uint32>(Direction::X);
    static constexpr uint32 idirY = static_cast<uint32>(Direction::Y);
    static constexpr uint32 idirZ = static_cast<uint32>(Direction::Z);

    static constexpr uint32 iBx = static_cast<uint32>(HybridQuantity::Scalar::Bx);
    static constexpr uint32 iBy = static_cast<uint32>(HybridQuantity::Scalar::By);
    static constexpr uint32 iBz = static_cast<uint32>(HybridQuantity::Scalar::Bz);

    static constexpr uint32 iEx = static_cast<uint32>(HybridQuantity::Scalar::Ex);
    static constexpr uint32 iEy = static_cast<uint32>(HybridQuantity::Scalar::Ey);
    static constexpr uint32 iEz = static_cast<uint32>(HybridQuantity::Scalar::Ez);

    static constexpr uint32 iJx = static_cast<uint32>(HybridQuantity::Scalar::Jx);
    static constexpr uint32 iJy = static_cast<uint32>(HybridQuantity::Scalar::Jy);
    static constexpr uint32 iJz = static_cast<uint32>(HybridQuantity::Scalar::Jz);

    static constexpr uint32 irho = static_cast<uint32>(HybridQuantity::Scalar::rho);

    static constexpr uint32 iVx = static_cast<uint32>(HybridQuantity::Scalar::Vx);
    static constexpr uint32 iVy = static_cast<uint32>(HybridQuantity::Scalar::Vy);
    static constexpr uint32 iVz = static_cast<uint32>(HybridQuantity::Scalar::Vz);

    static constexpr uint32 iP = static_cast<uint32>(HybridQuantity::Scalar::P);
};
} // namespace PHARE

#endif // GRIDLAYOUTDEFS_H
