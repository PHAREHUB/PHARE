#ifndef PHARE_CORE_HYBRID_HYBRID_QUANTITIES_HPP
#define PHARE_CORE_HYBRID_HYBRID_QUANTITIES_HPP

#include "core/def.hpp"
#include "core/utilities/enum_table.hpp"

#include <array>
#include <cstdint>
#include <tuple>
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
        count
    };
    enum class Vector : std::uint16_t { B = 0, E, J, V, count };
    enum class Tensor : std::uint16_t { M = 0, count };

    static constexpr auto vector_components = make_enum_table<Vector, std::array<Scalar, 3>>({
        {Vector::B, {Scalar::Bx, Scalar::By, Scalar::Bz}},
        {Vector::E, {Scalar::Ex, Scalar::Ey, Scalar::Ez}},
        {Vector::J, {Scalar::Jx, Scalar::Jy, Scalar::Jz}},
        {Vector::V, {Scalar::Vx, Scalar::Vy, Scalar::Vz}},
    });

    static constexpr auto all_primal_field = Scalar::rho;

    template<std::size_t rank, typename = std::enable_if_t<rank == 1 or rank == 2, void>>
    using TensorType = std::conditional_t<rank == 1, Vector, Tensor>;

    NO_DISCARD static constexpr auto B() { return componentsQuantities(Vector::B); }
    NO_DISCARD static constexpr auto E() { return componentsQuantities(Vector::E); }
    NO_DISCARD static constexpr auto J() { return componentsQuantities(Vector::J); }
    NO_DISCARD static constexpr auto V() { return componentsQuantities(Vector::V); }

    NO_DISCARD static constexpr std::array<Scalar, 3> componentsQuantities(Vector const qty)
    {
        auto const idx = static_cast<std::size_t>(qty);
        if (idx >= vector_components.size())
            throw std::runtime_error("Error - invalid Vector");
        return vector_components[idx];
    }

    NO_DISCARD static constexpr std::array<Scalar, 6> componentsQuantities(Tensor qty)
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
