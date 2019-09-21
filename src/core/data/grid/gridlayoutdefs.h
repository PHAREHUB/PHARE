#ifndef PHARE_CORE_GRID_GRIDLAYOUTDEFS_H
#define PHARE_CORE_GRID_GRIDLAYOUTDEFS_H

#include <cstddef>

#include "hybrid/hybrid_quantities.h"
#include "utilities/point/point.h"
#include "utilities/types.h"

namespace PHARE
{
namespace core
{
    template<typename T>
    class StrongType
    {
    public:
        constexpr StrongType(T val)
            : value_(val)
        {
        }
        constexpr StrongType()
            : value_(0)
        {
        }
        constexpr operator T() const noexcept { return value_; }
        constexpr inline const T& value() const { return value_; }

    private:
        T value_;
    };

    struct Direction : public StrongType<uint32_t>
    {
        struct X_t;
        static X_t const X;
        struct Y_t;
        static Y_t const Y;
        struct Z_t;
        static Z_t const Z;
        constexpr bool operator==(const Direction& that) { return this->value() == that.value(); }
        constexpr operator uint32_t() const noexcept { return value(); }
    };
    struct Direction::X_t : public Direction
    {
    };
    inline Direction::X_t constexpr const Direction::X{{0}};

    struct Direction::Y_t : public Direction
    {
    };
    inline Direction::Y_t constexpr const Direction::Y{{1}};

    struct Direction::Z_t : public Direction
    {
    };
    inline Direction::Z_t constexpr const Direction::Z{{2}};

    struct QtyCentering : public StrongType<uint32_t>
    {
        struct primal_t;
        static primal_t const primal;
        struct dual_t;
        static dual_t const dual;
        constexpr bool operator==(const QtyCentering& that)
        {
            return this->value() == that.value();
        }
        constexpr operator uint32_t() const noexcept { return value(); }
    };
    struct QtyCentering::primal_t : public QtyCentering
    {
    };
    inline QtyCentering::primal_t constexpr const QtyCentering::primal{{0}};

    struct QtyCentering::dual_t : public QtyCentering
    {
    };
    inline QtyCentering::dual_t constexpr const QtyCentering::dual{{1}};

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

        static constexpr uint32 idirX = Direction::X;
        static constexpr uint32 idirY = Direction::Y;
        static constexpr uint32 idirZ = Direction::Z;

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
} // namespace core
} // namespace PHARE

#endif // GRIDLAYOUTDEFS_H
