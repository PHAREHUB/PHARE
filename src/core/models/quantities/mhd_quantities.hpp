#ifndef PHARE_CORE_MHD_MHD_QUANTITIES_HPP
#define PHARE_CORE_MHD_MHD_QUANTITIES_HPP

#include "core/def.hpp"
#include "core/utilities/enum_table.hpp"

#include <array>
#include <cstdint>
#include <tuple>
#include <stdexcept>


namespace PHARE::core
{
class MHDQuantity
{
public:
    enum class Scalar : std::uint16_t {
        rho = 0, // density
        Vx,      // velocity components
        Vy,
        Vz,
        Bx,
        By,
        Bz,
        P, // pressure

        Etot,  // total energy
        rhoVx, // momentum components
        rhoVy,
        rhoVz,

        Ex, // electric field components
        Ey,
        Ez,
        Jx, // current density components
        Jy,
        Jz,

        ScalarFlux_x,
        ScalarFlux_y,
        ScalarFlux_z,
        VecFluxX_x,
        VecFluxY_x,
        VecFluxZ_x,
        VecFluxX_y,
        VecFluxY_y,
        VecFluxZ_y,
        VecFluxX_z,
        VecFluxY_z,
        VecFluxZ_z,

        ScalarAllPrimal,
        VecAllPrimalX,
        VecAllPrimalY,
        VecAllPrimalZ,

        count
    };
    enum class Vector : std::uint16_t {
        V = 0,
        B,
        rhoV,
        E,
        J,
        VecFlux_x,
        VecFlux_y,
        VecFlux_z,
        VecAllPrimal,
        count
    };
    enum class Tensor : std::uint16_t { count = 0 };

    static constexpr auto vector_components = make_enum_table<Vector, std::array<Scalar, 3>>({
        {Vector::V, {Scalar::Vx, Scalar::Vy, Scalar::Vz}},
        {Vector::B, {Scalar::Bx, Scalar::By, Scalar::Bz}},
        {Vector::rhoV, {Scalar::rhoVx, Scalar::rhoVy, Scalar::rhoVz}},
        {Vector::E, {Scalar::Ex, Scalar::Ey, Scalar::Ez}},
        {Vector::J, {Scalar::Jx, Scalar::Jy, Scalar::Jz}},
        {Vector::VecFlux_x, {Scalar::VecFluxX_x, Scalar::VecFluxY_x, Scalar::VecFluxZ_x}},
        {Vector::VecFlux_y, {Scalar::VecFluxX_y, Scalar::VecFluxY_y, Scalar::VecFluxZ_y}},
        {Vector::VecFlux_z, {Scalar::VecFluxX_z, Scalar::VecFluxY_z, Scalar::VecFluxZ_z}},
        {Vector::VecAllPrimal,
         {Scalar::VecAllPrimalX, Scalar::VecAllPrimalY, Scalar::VecAllPrimalZ}},
    });

    static constexpr auto all_primal_field = Scalar::ScalarAllPrimal;

    template<std::size_t rank, typename = std::enable_if_t<rank == 1 or rank == 2, void>>
    using TensorType = std::conditional_t<rank == 1, Vector, Tensor>;

    NO_DISCARD static constexpr auto V() { return componentsQuantities(Vector::V); }
    NO_DISCARD static constexpr auto B() { return componentsQuantities(Vector::B); }
    NO_DISCARD static constexpr auto rhoV() { return componentsQuantities(Vector::rhoV); }

    NO_DISCARD static constexpr auto E() { return componentsQuantities(Vector::E); }
    NO_DISCARD static constexpr auto J() { return componentsQuantities(Vector::J); }

    NO_DISCARD static constexpr auto VecFlux_x() { return componentsQuantities(Vector::VecFlux_x); }
    NO_DISCARD static constexpr auto VecFlux_y() { return componentsQuantities(Vector::VecFlux_y); }
    NO_DISCARD static constexpr auto VecFlux_z() { return componentsQuantities(Vector::VecFlux_z); }

    NO_DISCARD static constexpr auto VecAllPrimal()
    {
        return componentsQuantities(Vector::VecAllPrimal);
    }

    NO_DISCARD static constexpr std::array<Scalar, 3> componentsQuantities(Vector const qty)
    {
        auto const idx = static_cast<std::size_t>(qty);
        if (idx >= vector_components.size())
            throw std::runtime_error("Error - invalid Vector");
        return vector_components[idx];
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
