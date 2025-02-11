#ifndef CORE_NUMERICS_GODUNOV_GODUNOV_UTILS_HPP
#define CORE_NUMERICS_GODUNOV_GODUNOV_UTILS_HPP

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"
#include "core/utilities/index/index.hpp"

namespace PHARE::core
{

struct PerIndex
{
    PerIndex(double rho, double Vx, double Vy, double Vz, double Bx, double By, double Bz, double P)
        : rho{rho}
        , Vx{Vx}
        , Vy{Vy}
        , Vz{Vz}
        , Bx{Bx}
        , By{By}
        , Bz{Bz}
        , P{P}
    {
    }

    auto operator()() const { return std::make_tuple(rho, Vx, Vy, Vz, Bx, By, Bz, P); }

    void toConservative(auto const& gamma)
    {
        auto const [rhoVx, rhoVy, rhoVz] = vToRhoV(rho, Vx, Vy, Vz);
        double Etot                      = eosPToEtot(gamma, rho, Vx, Vy, Vz, Bx, By, Bz, P);

        Vx = rhoVx;
        Vy = rhoVy;
        Vz = rhoVz;
        P  = Etot;
    }

    double rho;
    double Vx;
    double Vy;
    double Vz;
    double Bx;
    double By;
    double Bz;
    double P;
};

struct PerIndexRef
{
    PerIndexRef(double& rho, double& rhoVx, double& rhoVy, double& rhoVz, double& Bx, double& By,
                double& Bz, double& Etot)
        : rho{rho}
        , rhoVx{rhoVx}
        , rhoVy{rhoVy}
        , rhoVz{rhoVz}
        , Bx{Bx}
        , By{By}
        , Bz{Bz}
        , Etot{Etot}
    {
    }

    PerIndexRef& operator=(const PerIndex& other)
    {
        rho   = other.rho;
        rhoVx = other.Vx;
        rhoVy = other.Vy;
        rhoVz = other.Vz;
        Bx    = other.Bx;
        By    = other.By;
        Bz    = other.Bz;
        Etot  = other.P;
        return *this;
    }

    double& rho;
    double& rhoVx;
    double& rhoVy;
    double& rhoVz;
    double& Bx;
    double& By;
    double& Bz;
    double& Etot;
};

template<typename Field>
struct AllFluxes
{
    AllFluxes(auto&& fluxes)
        : rho_fx{std::get<0>(fluxes)}
        , rhoVx_fx{std::get<1>(fluxes)(Component::X)}
        , rhoVy_fx{std::get<1>(fluxes)(Component::Y)}
        , rhoVz_fx{std::get<1>(fluxes)(Component::Z)}
        , Bx_fx{std::get<2>(fluxes)(Component::X)}
        , By_fx{std::get<2>(fluxes)(Component::Y)}
        , Bz_fx{std::get<2>(fluxes)(Component::Z)}
        , Etot_fx{std::get<3>(fluxes)}
        , rho_fy{std::get<4>(fluxes)}
        , rhoVx_fy{std::get<5>(fluxes)(Component::X)}
        , rhoVy_fy{std::get<5>(fluxes)(Component::Y)}
        , rhoVz_fy{std::get<5>(fluxes)(Component::Z)}
        , Bx_fy{std::get<6>(fluxes)(Component::X)}
        , By_fy{std::get<6>(fluxes)(Component::Y)}
        , Bz_fy{std::get<6>(fluxes)(Component::Z)}
        , Etot_fy{std::get<7>(fluxes)}
        , rho_fz{std::get<8>(fluxes)}
        , rhoVx_fz{std::get<9>(fluxes)(Component::X)}
        , rhoVy_fz{std::get<9>(fluxes)(Component::Y)}
        , rhoVz_fz{std::get<9>(fluxes)(Component::Z)}
        , Bx_fz{std::get<10>(fluxes)(Component::X)}
        , By_fz{std::get<10>(fluxes)(Component::Y)}
        , Bz_fz{std::get<10>(fluxes)(Component::Z)}
        , Etot_fz{std::get<11>(fluxes)}
    {
    }

    template<auto direction>
    PerIndexRef get_dir(MeshIndex<Field::dimension> index) const
    {
        if constexpr (direction == Direction::X)
            return PerIndexRef{rho_fx(index), rhoVx_fx(index), rhoVy_fx(index), rhoVz_fx(index),
                               Bx_fx(index),  By_fx(index),    Bz_fx(index),    Etot_fx(index)};
        else if constexpr (direction == Direction::Y)
            return PerIndexRef{rho_fy(index), rhoVx_fy(index), rhoVy_fy(index), rhoVz_fy(index),
                               Bx_fy(index),  By_fy(index),    Bz_fy(index),    Etot_fy(index)};
        else if constexpr (direction == Direction::Z)
            return PerIndexRef{rho_fz(index), rhoVx_fz(index), rhoVy_fz(index), rhoVz_fz(index),
                               Bx_fz(index),  By_fz(index),    Bz_fz(index),    Etot_fz(index)};
    }

    Field& rho_fx;
    Field& rhoVx_fx;
    Field& rhoVy_fx;
    Field& rhoVz_fx;
    Field& Bx_fx;
    Field& By_fx;
    Field& Bz_fx;
    Field& Etot_fx;

    Field& rho_fy;
    Field& rhoVx_fy;
    Field& rhoVy_fy;
    Field& rhoVz_fy;
    Field& Bx_fy;
    Field& By_fy;
    Field& Bz_fy;
    Field& Etot_fy;

    Field& rho_fz;
    Field& rhoVx_fz;
    Field& rhoVy_fz;
    Field& rhoVz_fz;
    Field& Bx_fz;
    Field& By_fz;
    Field& Bz_fz;
    Field& Etot_fz;
};

} // namespace PHARE::core

#endif
