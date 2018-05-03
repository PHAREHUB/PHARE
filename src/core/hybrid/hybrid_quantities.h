#ifndef PHARE_CORE_HYBRID_HYBRID_QUANTITIES_H
#define PHARE_CORE_HYBRID_HYBRID_QUANTITIES_H



#include <array>

namespace PHARE
{
class HybridQuantity
{
public:
    enum class Scalar { Bx, By, Bz, Ex, Ey, Ez, Jx, Jy, Jz, rho, Vx, Vy, Vz, P, count };

    enum class Vector { B, E, J, V };

    static constexpr std::array<Scalar, 3> componentsQuantities(Vector qty)
    {
        if (qty == Vector::B)
        {
            return {{Scalar::Bx, Scalar::By, Scalar::Bz}};
        }

        else if (qty == Vector::E)
        {
            return {{Scalar::Ex, Scalar::Ey, Scalar::Ez}};
        }

        else if (qty == Vector::J)
        {
            return {{Scalar::Jx, Scalar::Jy, Scalar::Jz}};
        }

        else if (qty == Vector::V)
        {
            return {{Scalar::Vx, Scalar::Vy, Scalar::Vz}};
        }
    }
};
} // namespace PHARE

#endif
