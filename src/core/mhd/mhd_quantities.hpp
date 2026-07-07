#ifndef PHARE_CORE_MHD_MHD_QUANTITIES_HPP
#define PHARE_CORE_MHD_MHD_QUANTITIES_HPP

#include "core/def.hpp"

#include <array>
#include <tuple>
#include <stdexcept>


namespace PHARE::core
{
class MHDQuantity
{
public:
    enum class Scalar {
        rho, // density
        Vx,  // velocity components
        Vy,
        Vz,
        B1x, // perturbation magnetic field components (the evolved field)
        B1y,
        B1z,
        Bx, // total magnetic field components (B = B0 + B1, for diagnostics)
        By,
        Bz,
        B0x, // background magnetic field components (static, analytic)
        B0y,
        B0z,
        P, // pressure

        Etot1, // total energy using B1 only (the conserved energy variable)
        Etot,  // total energy using total B (for diagnostics)
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

        divB, // cell-centered divergence of the total field B, for diagnostics

        ScalarAllPrimal,
        VecAllPrimalX,
        VecAllPrimalY,
        VecAllPrimalZ,

        count
    };
    enum class Vector { V, B1, B, B0, rhoV, E, J, VecFlux_x, VecFlux_y, VecFlux_z, VecAllPrimal };
    enum class Tensor { count };

    static constexpr auto all_primal_field = Scalar::ScalarAllPrimal;

    template<std::size_t rank, typename = std::enable_if_t<rank == 1 or rank == 2, void>>
    using TensorType = std::conditional_t<rank == 1, Vector, Tensor>;

    NO_DISCARD static constexpr auto V() { return componentsQuantities(Vector::V); }
    NO_DISCARD static constexpr auto B1() { return componentsQuantities(Vector::B1); }
    NO_DISCARD static constexpr auto B() { return componentsQuantities(Vector::B); }
    NO_DISCARD static constexpr auto B0() { return componentsQuantities(Vector::B0); }
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

    NO_DISCARD static constexpr std::array<Scalar, 3> componentsQuantities(Vector qty)
    {
        if (qty == Vector::V)
            return {{Scalar::Vx, Scalar::Vy, Scalar::Vz}};

        if (qty == Vector::B1)
            return {{Scalar::B1x, Scalar::B1y, Scalar::B1z}};

        if (qty == Vector::B)
            return {{Scalar::Bx, Scalar::By, Scalar::Bz}};

        if (qty == Vector::B0)
            return {{Scalar::B0x, Scalar::B0y, Scalar::B0z}};

        if (qty == Vector::rhoV)
            return {{Scalar::rhoVx, Scalar::rhoVy, Scalar::rhoVz}};


        if (qty == Vector::E)
            return {{Scalar::Ex, Scalar::Ey, Scalar::Ez}};

        if (qty == Vector::J)
            return {{Scalar::Jx, Scalar::Jy, Scalar::Jz}};


        if (qty == Vector::VecFlux_x)
            return {{Scalar::VecFluxX_x, Scalar::VecFluxY_x, Scalar::VecFluxZ_x}};

        if (qty == Vector::VecFlux_y)
            return {{Scalar::VecFluxX_y, Scalar::VecFluxY_y, Scalar::VecFluxZ_y}};

        if (qty == Vector::VecFlux_z)
            return {{Scalar::VecFluxX_z, Scalar::VecFluxY_z, Scalar::VecFluxZ_z}};

        if (qty == Vector::VecAllPrimal)
            return {{Scalar::VecAllPrimalX, Scalar::VecAllPrimalY, Scalar::VecAllPrimalZ}};

        throw std::runtime_error("Error - invalid Vector");
    }

    NO_DISCARD static constexpr auto B1_items()
    {
        auto const& [B1x, B1y, B1z] = B1();
        return std::make_tuple(std::make_pair("B1x", B1x), std::make_pair("B1y", B1y),
                               std::make_pair("B1z", B1z));
    }
    NO_DISCARD static constexpr auto B_items()
    {
        auto const& [Bx, By, Bz] = B();
        return std::make_tuple(std::make_pair("Bx", Bx), std::make_pair("By", By),
                               std::make_pair("Bz", Bz));
    }
    NO_DISCARD static constexpr auto B0_items()
    {
        auto const& [B0x, B0y, B0z] = B0();
        return std::make_tuple(std::make_pair("B0x", B0x), std::make_pair("B0y", B0y),
                               std::make_pair("B0z", B0z));
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
