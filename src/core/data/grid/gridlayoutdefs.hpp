#ifndef PHARE_CORE_GRID_GRIDLAYOUTDEFS_HPP
#define PHARE_CORE_GRID_GRIDLAYOUTDEFS_HPP

#include <cstddef>

#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"

namespace PHARE
{
namespace core
{
    enum class Direction { X, Y, Z };


    enum class QtyCentering { primal = 0, dual = 1 };


    template<std::size_t dim>
    struct WeightPoint
    {
        constexpr WeightPoint(Point<int, dim> point, double _coef)
            : indexes{std::move(point)}
            , coef{_coef}
        {
        }

        Point<int, dim> indexes;
        double coef;
    };


    // using LinearCombination = std::vector<WeightPoint>;

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

        static constexpr std::uint32_t idirX = static_cast<std::uint32_t>(Direction::X);
        static constexpr std::uint32_t idirY = static_cast<std::uint32_t>(Direction::Y);
        static constexpr std::uint32_t idirZ = static_cast<std::uint32_t>(Direction::Z);

        static constexpr std::uint32_t iBx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Bx);
        static constexpr std::uint32_t iBy = static_cast<std::uint32_t>(HybridQuantity::Scalar::By);
        static constexpr std::uint32_t iBz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Bz);

        static constexpr std::uint32_t iEx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Ex);
        static constexpr std::uint32_t iEy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Ey);
        static constexpr std::uint32_t iEz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Ez);

        static constexpr std::uint32_t iJx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Jx);
        static constexpr std::uint32_t iJy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Jy);
        static constexpr std::uint32_t iJz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Jz);

        static constexpr std::uint32_t irho
            = static_cast<std::uint32_t>(HybridQuantity::Scalar::rho);

        static constexpr std::uint32_t iVx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Vx);
        static constexpr std::uint32_t iVy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Vy);
        static constexpr std::uint32_t iVz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Vz);

        static constexpr std::uint32_t iMxx
            = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mxx);
        static constexpr std::uint32_t iMxy
            = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mxy);
        static constexpr std::uint32_t iMxz
            = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mxz);
        static constexpr std::uint32_t iMyy
            = static_cast<std::uint32_t>(HybridQuantity::Scalar::Myy);
        static constexpr std::uint32_t iMyz
            = static_cast<std::uint32_t>(HybridQuantity::Scalar::Myz);
        static constexpr std::uint32_t iMzz
            = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mzz);

        static constexpr std::uint32_t iP = static_cast<std::uint32_t>(HybridQuantity::Scalar::P);
    };
} // namespace core
} // namespace PHARE

#endif // GRIDLAYOUTDEFS_HPP
