#ifndef PHARE_CORE_HYBRID_HYBRID_QUANTITIES_H
#define PHARE_CORE_HYBRID_HYBRID_QUANTITIES_H



#include <array>

namespace PHARE
{
class HybridQuantity
{
public:
    enum class Quantity {
        B,
        Bx,
        By,
        Bz,
        E,
        Ex,
        Ey,
        Ez,
        J,
        Jx,
        Jy,
        Jz,
        rho,
        V,
        Vx,
        Vy,
        Vz,
        P,
        count
    };

    static constexpr std::array<Quantity, 3> componentsQuantities(Quantity qty)
    {
        if (qty == Quantity::B)
        {
            return {{Quantity::Bx, Quantity::By, Quantity::Bz}};
        }

        else if (qty == Quantity::E)
        {
            return {{Quantity::Ex, Quantity::Ey, Quantity::Ez}};
        }

        else if (qty == Quantity::J)
        {
            return {{Quantity::Jx, Quantity::Jy, Quantity::Jz}};
        }

        else if (qty == Quantity::V)
        {
            return {{Quantity::Vx, Quantity::Vy, Quantity::Vz}};
        }
    }
};
} // namespace PHARE

#endif
