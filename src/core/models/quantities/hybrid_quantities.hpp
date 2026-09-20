#ifndef PHARE_CORE_HYBRID_HYBRID_QUANTITIES_HPP
#define PHARE_CORE_HYBRID_HYBRID_QUANTITIES_HPP

#include "core/def.hpp"


#include <array>
#include <cstdint>
#include <stdexcept>


namespace PHARE::core
{
class HybridQuantity
{
public:
    enum class Scalar : std::uint16_t {
        Bx = 0, // magnetic field components
        By,
        Bz,
        Ex, // electric field components
        Ey,
        Ez,
        Jx, // current density components
        Jy,
        Jz,
        rho, // charge density
        Vx,  // bulk velocity components
        Vy,
        Vz,
        P,   // pressure
        Mxx, // momentum tensor components
        Mxy,
        Mxz,
        Myy,
        Myz,
        Mzz,
        count,
        INVALID
    };

    enum class Vector { B, E, J, V, F, INVALID };
    enum class Tensor { M, count, INVALID };

    static constexpr auto all_primal_field = Scalar::rho;

    template<std::size_t rank, typename = std::enable_if_t<rank == 1 or rank == 2, void>>
    using TensorType = std::conditional_t<rank == 1, Vector, Tensor>;


    NO_DISCARD static constexpr auto B() { return componentsQuantities(Vector::B); }
    NO_DISCARD static constexpr auto E() { return componentsQuantities(Vector::E); }
    NO_DISCARD static constexpr auto J() { return componentsQuantities(Vector::J); }
    NO_DISCARD static constexpr auto V() { return componentsQuantities(Vector::V); }

    NO_DISCARD static constexpr std::array<Scalar, 3> componentsQuantities(Vector qty)
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

    NO_DISCARD static constexpr std::array<Scalar, 6> componentsQuantities(Tensor /*qty*/)
    {
        // no condition, for now there's only then momentum tensor M
        return {{Scalar::Mxx, Scalar::Mxy, Scalar::Mxz, Scalar::Myy, Scalar::Myz, Scalar::Mzz}};
    }
};

auto inline componentsQuantities(HybridQuantity::Vector const qty)
{
    return HybridQuantity::componentsQuantities(qty);
}

auto inline componentsQuantities(HybridQuantity::Tensor const qty)
{
    return HybridQuantity::componentsQuantities(qty);
}

} // namespace PHARE::core

#endif
