#ifndef PHARE_CORE_PIC_PIC_QUANTITIES_HPP
#define PHARE_CORE_PIC_PIC_QUANTITIES_HPP

#include "core/def.hpp"

#include <array>
#include <tuple>
#include <stdexcept>


namespace PHARE::core
{
class PICQuantity
{
public:
    enum class Scalar {
        Bx, // magnetic field components
        By,
        Bz,
        Ex, // electric field components
        Ey,
        Ez,
        Jx, // current density components
        Jy,
        Jz,
        rho, // charge density
        Vx,  // bulk velocity components for ions
        Vy,
        Vz,
        Mxx, // momentum tensor components
        Mxy,
        Mxz,
        Myy,
        Myz,
        Mzz,
        Vex,  // bulk velocity components for electrons
        Vey,
        Vez,
        count
    };
    enum class Vector { B, E, J, V, Ve };
    enum class Tensor { M, count };

    template<std::size_t rank, typename = std::enable_if_t<rank == 1 or rank == 2, void>>
    using TensorType = std::conditional_t<rank == 1, Vector, Tensor>;

    NO_DISCARD static constexpr auto B() { return componentsQuantities(Vector::B); }
    NO_DISCARD static constexpr auto E() { return componentsQuantities(Vector::E); }
    NO_DISCARD static constexpr auto J() { return componentsQuantities(Vector::J); }
    NO_DISCARD static constexpr auto V() { return componentsQuantities(Vector::V); }
    NO_DISCARD static constexpr auto Ve() { return componentsQuantities(Vector::Ve); }

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

        if (qty == Vector::Ve)
            return {{Scalar::Vex, Scalar::Vey, Scalar::Vez}};

        throw std::runtime_error("Error - invalid Vector");
    }
    static constexpr std::array<Scalar, 6> componentsQuantities(Tensor qty)
    {
        // no condition, for now there's only then momentum tensor M
        return {{Scalar::Mxx, Scalar::Mxy, Scalar::Mxz, Scalar::Myy, Scalar::Myz, Scalar::Mzz}};
    }

    NO_DISCARD static constexpr auto B_items()
    {
        auto const& [Bx, By, Bz] = B();
        return std::make_tuple(std::make_pair("Bx", Bx), std::make_pair("By", By),
                               std::make_pair("Bz", Bz));
    }
    NO_DISCARD static constexpr auto E_items()
    {
        auto const& [Ex, Ey, Ez] = E();
        return std::make_tuple(std::make_pair("Ex", Ex), std::make_pair("Ey", Ey),
                               std::make_pair("Ez", Ez));
    }
};

} // namespace PHARE::core

#endif
