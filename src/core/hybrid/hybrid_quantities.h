#ifndef PHARE_CORE_HYBRID_HYBRID_QUANTITIES_H
#define PHARE_CORE_HYBRID_HYBRID_QUANTITIES_H


#include <stdexcept>
#include <array>

namespace PHARE
{
namespace core
{
    class HybridQuantity
    {
    public:
        enum class Scalar {
            Bx,
            By,
            Bz,
            Ex,
            Ey,
            Ez,
            Jx,
            Jy,
            Jz,
            rho,
            Vx,
            Vy,
            Vz,
            P,
            count,
            INVALID
        };

        enum class Vector { B, E, J, V, INVALID };

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
            else
            {
                throw std::runtime_error("Error - invalid Vector");
            }
        }
    };
} // namespace core
} // namespace PHARE

#endif
