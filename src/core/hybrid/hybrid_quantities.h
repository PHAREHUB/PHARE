#ifndef PHARE_CORE_HYBRID_HYBRID_QUANTITIES_H
#define PHARE_CORE_HYBRID_HYBRID_QUANTITIES_H


#include <array>
#include <tuple>
#include <stdexcept>


namespace PHARE::core
{
class HybridQuantity
{
public:
    enum class Scalar { Bx, By, Bz, Ex, Ey, Ez, Jx, Jy, Jz, rho, Vx, Vy, Vz, P, count };
    enum class Vector { B, E, J, V };

    static constexpr auto B() { return componentsQuantities(Vector::B); }
    static constexpr auto E() { return componentsQuantities(Vector::E); }
    static constexpr auto J() { return componentsQuantities(Vector::J); }
    static constexpr auto V() { return componentsQuantities(Vector::V); }

    static constexpr std::array<Scalar, 3> componentsQuantities(Vector qty)
    {
        if (qty == Vector::B)
            return {{Scalar::Bx, Scalar::By, Scalar::Bz}};

        if (qty == Vector::E)
            return {{Scalar::Ex, Scalar::Ey, Scalar::Ez}};

        if (qty == Vector::J)
            return {{Scalar::Jx, Scalar::Jy, Scalar::Jz}};

        if (qty == Vector::V)
            return {{Scalar::Vx, Scalar::Vy, Scalar::Vz}};

        throw std::runtime_error("Error - invalid Vector");
    }

    static constexpr auto B_items()
    {
        auto const& [Bx, By, Bz] = B();
        return std::make_tuple(std::make_pair("Bx", Bx), std::make_pair("By", By),
                               std::make_pair("Bz", Bz));
    }
    static constexpr auto E_items()
    {
        auto const& [Ex, Ey, Ez] = E();
        return std::make_tuple(std::make_pair("Ex", Ex), std::make_pair("Ey", Ey),
                               std::make_pair("Ez", Ez));
    }
};

} // namespace PHARE::core

namespace std
{
auto inline to_string(PHARE::core::HybridQuantity::Scalar v)
{
    return std::to_string(static_cast<int>(v));
}
}

#endif
