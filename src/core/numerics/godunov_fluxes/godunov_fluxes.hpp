#ifndef PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP
#define PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include "initializer/data_provider.hpp"
#include <cstddef>
#include <tuple>

namespace PHARE::core
{
template<typename GridLayout>
class GodunovFluxes : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    GodunovFluxes(PHARE::initializer::PHAREDict const& dict)
        : gamma_{dict["heat_capacity_ratio"].template to<double>()}
        , eta_{dict["resistivity"].template to<double>()}
        , nu_{dict["hyper_resistivity"].template to<double>()}

        , terms_{dict["terms"].template to<std::string>()}
        , reconstruction_{dict["reconstruction"].template to<std::string>()}
        , limiter_{dict["limiter"].template to<std::string>()}
        , riemann_{dict["riemann"].template to<std::string>()}
    {
    }

    template<typename Field, typename VecField, typename... Fluxes>
    void operator()(Field const& rho, VecField const& V, VecField const& B, Field const& P,
                    VecField const& J, Fluxes&... fluxes) const
    {
        if (!this->hasLayout())
            throw std::runtime_error("Error - GodunovFluxes - GridLayout not set, cannot proceed "
                                     "to reconstruction");

        if constexpr (dimension == 1)
        {
            auto&& flux_tuple = std::forward_as_tuple(fluxes...);

            auto& rho_x  = std::get<0>(flux_tuple);
            auto& rhoV_x = std::get<1>(flux_tuple);
            auto& B_x    = std::get<2>(flux_tuple);
            auto& Etot_x = std::get<3>(flux_tuple);

            layout_->evalOnBox(rho_x, [&](auto&... args) mutable {
                this->template godunov_fluxes_<Direction::X>(rho, V, B, P, J, rho_x, rhoV_x, B_x,
                                                             Etot_x, {args...});
            });
        }
        if constexpr (dimension == 2)
        {
            auto&& flux_tuple = std::forward_as_tuple(fluxes...);

            auto& rho_x  = std::get<0>(flux_tuple);
            auto& rhoV_x = std::get<1>(flux_tuple);
            auto& B_x    = std::get<2>(flux_tuple);
            auto& Etot_x = std::get<3>(flux_tuple);

            auto& rho_y  = std::get<4>(flux_tuple);
            auto& rhoV_y = std::get<5>(flux_tuple);
            auto& B_y    = std::get<6>(flux_tuple);
            auto& Etot_y = std::get<7>(flux_tuple);

            layout_->evalOnBox(rho_x, [&](auto&... args) mutable {
                this->template godunov_fluxes_<Direction::X>(rho, V, B, P, J, rho_x, rhoV_x, B_x,
                                                             Etot_x, {args...});
            });

            layout_->evalOnBox(rho_y, [&](auto&... args) mutable {
                this->template godunov_fluxes_<Direction::Y>(rho, V, B, P, J, rho_y, rhoV_y, B_y,
                                                             Etot_y, {args...});
            });
        }
        if constexpr (dimension == 3)
        {
            auto&& flux_tuple = std::forward_as_tuple(fluxes...);

            auto& rho_x  = std::get<0>(flux_tuple);
            auto& rhoV_x = std::get<1>(flux_tuple);
            auto& B_x    = std::get<2>(flux_tuple);
            auto& Etot_x = std::get<3>(flux_tuple);

            auto& rho_y  = std::get<4>(flux_tuple);
            auto& rhoV_y = std::get<5>(flux_tuple);
            auto& B_y    = std::get<6>(flux_tuple);
            auto& Etot_y = std::get<7>(flux_tuple);

            auto& rho_z  = std::get<8>(flux_tuple);
            auto& rhoV_z = std::get<9>(flux_tuple);
            auto& B_z    = std::get<10>(flux_tuple);
            auto& Etot_z = std::get<11>(flux_tuple);

            layout_->evalOnBox(rho_x, [&](auto&... args) mutable {
                this->template godunov_fluxes_<Direction::X>(rho, V, B, P, J, rho_x, rhoV_x, B_x,
                                                             Etot_x, {args...});
            });

            layout_->evalOnBox(rho_y, [&](auto&... args) mutable {
                this->template godunov_fluxes_<Direction::Y>(rho, V, B, P, J, rho_y, rhoV_y, B_y,
                                                             Etot_y, {args...});
            });

            layout_->evalOnBox(rho_z, [&](auto&... args) mutable {
                this->template godunov_fluxes_<Direction::Z>(rho, V, B, P, J, rho_z, rhoV_z, B_z,
                                                             Etot_z, {args...});
            });
        }
    }

private:
    double const gamma_;
    double const eta_;
    double const nu_;

    std::string const terms_;
    std::string const reconstruction_;
    std::string const limiter_;
    std::string const riemann_;

    template<auto direction, typename Field, typename VecField>
    void godunov_fluxes_(Field const& rho, VecField const& V, VecField const& B, Field const& P,
                         VecField const& J, Field& rho_flux, VecField& rhoV_flux, VecField& B_flux,
                         Field& Etot_flux, MeshIndex<Field::dimension> index) const
    {
        auto const& [Vx, Vy, Vz] = V();
        auto const& [Bx, By, Bz] = B();
        auto const& [Jx, Jy, Jz] = J();

        auto [BxL, ByL, BzL, BxR, ByR, BzR] = center_reconstruct_<direction>(
            B, index, GridLayout::faceXToCellCenter(), GridLayout::faceYToCellCenter(),
            GridLayout::faceZToCellCenter());

        auto [JxL, JyL, JzL, JxR, JyR, JzR] = center_reconstruct_<direction>(
            J, index, GridLayout::edgeXToCellCenter(), GridLayout::edgeYToCellCenter(),
            GridLayout::edgeZToCellCenter());

        auto const& Quantities = std::forward_as_tuple(rho, Vx, Vy, Vz, P);

        // Left and Right state reconstructions
        auto [rhoL, VxL, VyL, VzL, PL] = std::apply(
            [&](auto const&... fields) {
                return std::make_tuple(reconstruct_uL_<direction>(fields, index)...);
            },
            Quantities);

        auto [rhoR, VxR, VyR, VzR, PR] = std::apply(
            [&](auto const&... fields) {
                return std::make_tuple(reconstruct_uR_<direction>(fields, index)...);
            },
            Quantities);

        // Compute ideal flux vector for Left and Right states
        auto [F_rhoL, F_rhoVxL, F_rhoVyL, F_rhoVzL, F_BxL, F_ByL, F_BzL, F_EtotL]
            = ideal_flux_vector_<direction>(rhoL, VxL, VyL, VzL, BxL, ByL, BzL, PL);

        auto [F_rhoR, F_rhoVxR, F_rhoVyR, F_rhoVzR, F_BxR, F_ByR, F_BzR, F_EtotR]
            = ideal_flux_vector_<direction>(rhoR, VxR, VyR, VzR, BxR, ByR, BzR, PR);

        if (terms_ == "hall")
        {
            hall_contribution_<direction>(rhoL, BxL, ByL, BzL, JxL, JyL, JzL, F_BxL, F_ByL, F_BzL,
                                          F_EtotL);

            hall_contribution_<direction>(rhoR, BxR, ByR, BzR, JxR, JyR, JzR, F_BxR, F_ByR, F_BzR,
                                          F_EtotR);

            /*resistive_contributions_<direction>(eta_, JxL, JyL, JzL, BxL, ByL, BzL, F_BxL,
             * F_ByL,*/
            /*                                    F_BzL, F_EtotL);*/
            /**/
            /*resistive_contributions_<direction>(eta_, JxR, JyR, JzR, BxR, ByR, BzR, F_BxR,
             * F_ByR,*/
            /*                                    F_BzR, F_EtotR);*/
            /**/
            /*auto [LaplJxL, LaplJxR] = reconstructed_laplacian(Jx, index);*/
            /*auto [LaplJyL, LaplJyR] = reconstructed_laplacian(Jy, index);*/
            /*auto [LaplJzL, LaplJzR] = reconstructed_laplacian(Jz, index);*/
            /**/
            /*resistive_contributions_<direction>(nu_, LaplJxL, LaplJyL, LaplJzL, BxL, ByL, BzL,*/
            /*                                    F_BxL, F_ByL, F_BzL, F_EtotL);*/
            /**/
            /*resistive_contributions_<direction>(nu_, LaplJxR, LaplJyR, LaplJzR, BxR, ByR, BzR,*/
            /*                                    F_BxR, F_ByR, F_BzR, F_EtotR);*/
        }

        auto uL = std::forward_as_tuple(rhoL, VxL, VyL, VzL, BxL, ByL, BzL, PL);
        auto uR = std::forward_as_tuple(rhoR, VxR, VyR, VzR, BxR, ByR, BzR, PR);
        auto fL = std::forward_as_tuple(F_rhoL, F_rhoVxL, F_rhoVyL, F_rhoVzL, F_BxL, F_ByL, F_BzL,
                                        F_EtotL);
        auto fR = std::forward_as_tuple(F_rhoR, F_rhoVxR, F_rhoVyR, F_rhoVzR, F_BxR, F_ByR, F_BzR,
                                        F_EtotR);

        auto [rho_, rhoVx_, rhoVy_, rhoVz_, Bx_, By_, Bz_, Etot_]
            = riemann_fluxes_<direction>(uL, uR, fL, fR);

        rho_flux(index)                = rho_;
        rhoV_flux(Component::X)(index) = rhoVx_;
        rhoV_flux(Component::Y)(index) = rhoVy_;
        rhoV_flux(Component::Z)(index) = rhoVz_;
        B_flux(Component::X)(index)    = Bx_;
        B_flux(Component::Y)(index)    = By_;
        B_flux(Component::Z)(index)    = Bz_;
        Etot_flux(index)               = Etot_;
    }

    template<auto direction>
    auto riemann_fluxes_(auto const& uL, auto const& uR, auto const& fL, auto const& fR) const
    {
        auto const& [rhoL, VxL, VyL, VzL, BxL, ByL, BzL, PL]                             = uL;
        auto const& [rhoR, VxR, VyR, VzR, BxR, ByR, BzR, PR]                             = uR;
        auto const& [F_rhoL, F_rhoVxL, F_rhoVyL, F_rhoVzL, F_BxL, F_ByL, F_BzL, F_EtotL] = fL;
        auto const& [F_rhoR, F_rhoVxR, F_rhoVyR, F_rhoVzR, F_BxR, F_ByR, F_BzR, F_EtotR] = fR;

        // Convert to conserved variables
        auto rhoVxL = rhoL * VxL;
        auto rhoVyL = rhoL * VyL;
        auto rhoVzL = rhoL * VzL;
        auto EtotL  = PL / (gamma_ - 1) + 0.5 * rhoL * (VxL * VxL + VyL * VyL + VzL * VzL)
                     + 0.5 * (BxL * BxL + ByL * ByL + BzL * BzL);

        auto rhoVxR = rhoR * VxR;
        auto rhoVyR = rhoR * VyR;
        auto rhoVzR = rhoR * VzR;
        auto EtotR  = PR / (gamma_ - 1) + 0.5 * rhoR * (VxR * VxR + VyR * VyR + VzR * VzR)
                     + 0.5 * (BxR * BxR + ByR * ByR + BzR * BzR);

        auto uL_ = std::forward_as_tuple(rhoL, rhoVxL, rhoVyL, rhoVzL);
        auto uR_ = std::forward_as_tuple(rhoR, rhoVxR, rhoVyR, rhoVzR);
        auto fL_ = std::forward_as_tuple(F_rhoL, F_rhoVxL, F_rhoVyL, F_rhoVzL);
        auto fR_ = std::forward_as_tuple(F_rhoR, F_rhoVxR, F_rhoVyR, F_rhoVzR);

        auto ubL = std::forward_as_tuple(BxL, ByL, BzL, EtotL);
        auto ubR = std::forward_as_tuple(BxR, ByR, BzR, EtotR);
        auto fbL = std::forward_as_tuple(F_BxL, F_ByL, F_BzL, F_EtotL);
        auto fbR = std::forward_as_tuple(F_BxR, F_ByR, F_BzR, F_EtotR);

        return riemann_steps_<direction>(uL_, ubL, uR_, ubR, fL_, fbL, fR_, fbR);
    }

    template<auto direction>
    auto riemann_steps_(auto const& uL, auto const& ubL, auto const& uR, auto const& ubR,
                        auto const& fL, auto const& fbL, auto const& fR, auto const& fbR) const
    {
        auto [speeds, speedsb]
            = riemann_speeds_<direction>(std::tuple_cat(uL, ubL), std::tuple_cat(uR, ubR));

        auto [rho, rhoVx, rhoVy, rhoVz] = std::apply(
            [&](auto... args) { return riemann_solver_(uL, uR, fL, fR, args...); }, speeds);

        auto [Bx, By, Bz, Etot] = std::apply(
            [&](auto... args) { return riemann_solver_(ubL, ubR, fbL, fbR, args...); }, speedsb);

        return std::make_tuple(rho, rhoVx, rhoVy, rhoVz, Bx, By, Bz, Etot);
    }

    template<auto direction>
    auto riemann_speeds_(auto const& uL, auto const& uR) const
    {
        if (riemann_ == "rusanov")
            return rusanov_speeds_<direction>(uL, uR);
        else
            throw std::runtime_error("Error - Riemann Solver - Unknown Riemann solver");
    }

    template<auto direction>
    auto rusanov_speeds_(auto const& uL, auto const& uR) const
    {
        auto const& [rhoL, VxL, VyL, VzL, BxL, ByL, BzL, PL] = uL;
        auto const& [rhoR, VxR, VyR, VzR, BxR, ByR, BzR, PR] = uR;
        auto BdotBL                                          = BxL * BxL + ByL * ByL + BzL * BzL;
        auto BdotBR                                          = BxR * BxR + ByR * ByR + BzR * BzR;

        auto compute_speeds = [&](auto rhoL, auto rhoR, auto PL, auto PR, auto BdotBL, auto BdotBR,
                                  auto VcompL, auto VcompR, auto BcompL, auto BcompR) {
            auto cfastL = compute_fast_magnetosonic_(rhoL, BcompL, BdotBL, PL);
            auto cfastR = compute_fast_magnetosonic_(rhoR, BcompR, BdotBR, PR);
            auto S      = std::max(std::abs(VcompL) + cfastL, std::abs(VcompR) + cfastR);
            auto Sb     = S;

            if (terms_ == "hall")
            {
                auto cwL = compute_whistler_(layout_->inverseMeshSize(direction), rhoL, BdotBL);
                auto cwR = compute_whistler_(layout_->inverseMeshSize(direction), rhoR, BdotBR);
                Sb = std::max(std::abs(VcompL) + cfastL + cwL, std::abs(VcompR) + cfastR + cwR);
            }

            return std::make_tuple(std::make_tuple(S), std::make_tuple(Sb));
        };

        if constexpr (direction == Direction::X)
            return compute_speeds(rhoL, rhoR, PL, PR, BdotBL, BdotBR, VxL, VxR, BxL, BxR);
        else if constexpr (direction == Direction::Y)
            return compute_speeds(rhoL, rhoR, PL, PR, BdotBL, BdotBR, VyL, VyR, ByL, ByR);
        else if constexpr (direction == Direction::Z)
            return compute_speeds(rhoL, rhoR, PL, PR, BdotBL, BdotBR, VzL, VzR, BzL, BzR);
    }

    template<typename... Speeds>
    auto riemann_solver_(auto const& uL, auto const& uR, auto const& fL, auto const& fR,
                         Speeds... S) const
    {
        if (riemann_ == "rusanov")
            return rusanov_(uL, uR, fL, fR, S...);
        else
            throw std::runtime_error("Error - Riemann Solver - Unknown Riemann solver");
    }

    auto rusanov_(auto const& uL, auto const& uR, auto const& fL, auto const& fR,
                  auto const S) const
    {
        // to be used 2 times in hall mhd (the second time for B with whisler contribution).

        auto constexpr N_elements = std::tuple_size_v<std::decay_t<decltype(uL)>>;

        return for_N<N_elements, for_N_R_mode::make_array>([&](auto i) {
            return (std::get<i>(fL) + std::get<i>(fR)) * 0.5
                   - S * (std::get<i>(uR) - std::get<i>(uL)) * 0.5;
        });
    }

    auto compute_fast_magnetosonic_(auto const& rho, auto const& B, auto const& BdotB,
                                    auto const& P) const
    {
        auto Sound     = std::sqrt((gamma_ * P) / rho);
        auto AlfvenDir = std::sqrt(B * B / rho); // directionnal alfven
        auto Alfven    = std::sqrt(BdotB / rho);

        auto c02    = Sound * Sound;
        auto cA2    = Alfven * Alfven;
        auto cAdir2 = AlfvenDir * AlfvenDir;

        return std::sqrt((c02 + cA2) * 0.5
                         + std::sqrt((c02 + cA2) * (c02 + cA2) - 4.0 * c02 * cAdir2) * 0.5);
    }

    auto compute_whistler_(auto const& invMeshSize, auto const& rho, auto const& BdotB) const
    {
        auto vw = std::sqrt(1 + 0.25 * invMeshSize * invMeshSize) + 0.5 * invMeshSize;
        return std::sqrt(BdotB) * vw / rho;
    }

    template<auto direction>
    auto ideal_flux_vector_(auto const& rho, auto const& Vx, auto const& Vy, auto const& Vz,
                            auto const& Bx, auto const& By, auto const& Bz, auto const& P) const
    {
        auto GeneralisedPressure = P + 0.5 * (Bx * Bx + By * By + Bz * Bz);
        auto TotalEnergy         = P / (gamma_ - 1) + 0.5 * rho * (Vx * Vx + Vy * Vy + Vz * Vz)
                           + 0.5 * (Bx * Bx + By * By + Bz * Bz);
        if constexpr (direction == Direction::X)
        {
            auto F_rho   = rho * Vx;
            auto F_rhoVx = rho * Vx * Vx + GeneralisedPressure - Bx * Bx;
            auto F_rhoVy = rho * Vx * Vy - Bx * By;
            auto F_rhoVz = rho * Vx * Vz - Bx * Bz;
            auto F_Bx    = 0.0;
            auto F_By    = By * Vx - Vy * Bx;
            auto F_Bz    = Bz * Vx - Vz * Bx;
            auto F_Etot
                = (TotalEnergy + GeneralisedPressure) * Vx - Bx * (Vx * Bx + Vy * By + Vz * Bz);

            return std::make_tuple(F_rho, F_rhoVx, F_rhoVy, F_rhoVz, F_Bx, F_By, F_Bz, F_Etot);
        }
        if constexpr (direction == Direction::Y)
        {
            auto F_rho   = rho * Vy;
            auto F_rhoVx = rho * Vy * Vx - By * Bx;
            auto F_rhoVy = rho * Vy * Vy + GeneralisedPressure - By * By;
            auto F_rhoVz = rho * Vy * Vz - By * Bz;
            auto F_Bx    = Bx * Vy - Vx * By;
            auto F_By    = 0.0;
            auto F_Bz    = Bz * Vy - Vz * By;
            auto F_Etot
                = (TotalEnergy + GeneralisedPressure) * Vy - By * (Vx * Bx + Vy * By + Vz * Bz);

            return std::make_tuple(F_rho, F_rhoVx, F_rhoVy, F_rhoVz, F_Bx, F_By, F_Bz, F_Etot);
        }
        if constexpr (direction == Direction::Z)
        {
            auto F_rho   = rho * Vz;
            auto F_rhoVx = rho * Vz * Vx - Bz * Bx;
            auto F_rhoVy = rho * Vz * Vy - Bz * By;
            auto F_rhoVz = rho * Vz * Vz + GeneralisedPressure - Bz * Bz;
            auto F_Bx    = Bx * Vz - Vx * Bz;
            auto F_By    = By * Vz - Vy * Bz;
            auto F_Bz    = 0.0;
            auto F_Etot
                = (TotalEnergy + GeneralisedPressure) * Vz - Bz * (Vx * Bx + Vy * By + Vz * Bz);

            return std::make_tuple(F_rho, F_rhoVx, F_rhoVy, F_rhoVz, F_Bx, F_By, F_Bz, F_Etot);
        }
    }

    template<auto direction>
    void hall_contribution_(auto const& rho, auto const& Bx, auto const& By, auto const& Bz,
                            auto const& Jx, auto const& Jy, auto const& Jz, auto& F_Bx, auto& F_By,
                            auto& F_Bz, auto& F_Etot) const
    {
        auto invRho = 1.0 / rho;

        auto JxB_x = Jy * Bz - Jz * By;
        auto JxB_y = Jz * Bx - Jx * Bz;
        auto JxB_z = Jx * By - Jy * Bx;

        auto BdotJ = Bx * Jx + By * Jy + Bz * Jz;
        auto BdotB = Bx * Bx + By * By + Bz * Bz;

        if constexpr (direction == Direction::X)
        {
            F_By += -JxB_z * invRho;
            F_Bz += JxB_y * invRho;
            F_Etot += (BdotJ * Bx - BdotB * Jx) * invRho;
        }
        if constexpr (direction == Direction::Y)
        {
            F_Bx += JxB_z * invRho;
            F_Bz += -JxB_x * invRho;
            F_Etot += (BdotJ * By - BdotB * Jy) * invRho;
        }
        if constexpr (direction == Direction::Z)
        {
            F_Bx += -JxB_y * invRho;
            F_By += JxB_x * invRho;
            F_Etot += (BdotJ * Bz - BdotB * Jz) * invRho;
        }
    }

    template<auto direction>
    void resistive_contributions_(auto const& pc, auto const& Jx, auto const& Jy, auto const& Jz,
                                  auto const& Bx, auto const& By, auto const& Bz, auto& F_Bx,
                                  auto& F_By, auto& F_Bz, auto& F_Etot) const
    // Can be used for both resistivity with J and eta and hyper resistivity with laplJ and nu
    {
        if constexpr (direction == Direction::X)
        {
            F_By += -Jz * pc;
            F_Bz += Jy * pc;
            F_Etot += (Jy * Bz - Jz * By) * pc;
        }
        if constexpr (direction == Direction::Y)
        {
            F_Bx += Jz * pc;
            F_Bz += -Jx * pc;
            F_Etot += (Jz * Bx - Jx * Bz) * pc;
        }
        if constexpr (direction == Direction::Y)
        {
            F_Bx += -Jy * pc;
            F_By += Jx * pc;
            F_Etot += (Jx * By - Jy * Bx) * pc;
        }
    }

    template<typename Field>
    auto reconstructed_laplacian(Field const& J, MeshIndex<Field::dimension> index) const
    {
        auto d2
            = [&](Direction dir, auto const& prevValue, auto const& Value, auto const& nextValue) {
                  return (layout_->inverseMeshSize(dir)) * (layout_->inverseMeshSize(dir))
                         * (prevValue - 2.0 * Value + nextValue);
              };

        auto JL = reconstruct_uL_<Direction::X>(J, index);
        auto JR = reconstruct_uR_<Direction::X>(J, index);

        if constexpr (dimension == 1)
        {
            MeshIndex<1> prevX = index;
            prevX[0] -= 1;
            MeshIndex<1> nextX = index;
            nextX[0] += 1;

            auto JL_X_1 = reconstruct_uL_<Direction::X>(J, prevX);
            auto JR_X_1 = reconstruct_uR_<Direction::X>(J, prevX);
            auto JL_X1  = reconstruct_uL_<Direction::X>(J, nextX);
            auto JR_X1  = reconstruct_uR_<Direction::X>(J, nextX);

            auto LaplJL = d2(Direction::X, JL_X_1, JL, JL_X1);
            auto LaplJR = d2(Direction::X, JR_X_1, JR, JR_X1);

            return std::make_tuple(LaplJL, LaplJR);
        }
        if constexpr (dimension == 2)
        {
            MeshIndex<2> prevX = index;
            prevX[0] -= 1;
            MeshIndex<2> nextX = index;
            nextX[0] += 1;

            MeshIndex<2> prevY = index;
            prevY[1] -= 1;
            MeshIndex<2> nextY = index;
            nextY[1] += 1;

            auto JL_X_1 = reconstruct_uL_<Direction::X>(J, prevX);
            auto JR_X_1 = reconstruct_uR_<Direction::X>(J, prevX);
            auto JL_X1  = reconstruct_uL_<Direction::X>(J, nextX);
            auto JR_X1  = reconstruct_uR_<Direction::X>(J, nextX);

            auto JL_Y_1 = reconstruct_uL_<Direction::Y>(J, prevY);
            auto JR_Y_1 = reconstruct_uR_<Direction::Y>(J, prevY);
            auto JL_Y1  = reconstruct_uL_<Direction::Y>(J, nextY);
            auto JR_Y1  = reconstruct_uR_<Direction::Y>(J, nextY);

            auto LaplJL = d2(Direction::X, JL_X_1, JL, JL_X1) + d2(Direction::Y, JL_Y_1, JL, JL_Y1);
            auto LaplJR = d2(Direction::X, JR_X_1, JR, JR_X1) + d2(Direction::Y, JR_Y_1, JR, JR_Y1);

            return std::make_tuple(LaplJL, LaplJR);
        }
        if constexpr (dimension == 3)
        {
            MeshIndex<3> prevX = index;
            prevX[0] -= 1;
            MeshIndex<3> nextX = index;
            nextX[0] += 1;

            MeshIndex<3> prevY = index;
            prevY[1] -= 1;
            MeshIndex<3> nextY = index;
            nextY[1] += 1;

            MeshIndex<3> prevZ = index;
            prevZ[2] -= 1;
            MeshIndex<3> nextZ = index;
            nextZ[2] += 1;

            auto JL_X_1 = reconstruct_uL_<Direction::X>(J, prevX);
            auto JR_X_1 = reconstruct_uR_<Direction::X>(J, prevX);
            auto JL_X1  = reconstruct_uL_<Direction::X>(J, nextX);
            auto JR_X1  = reconstruct_uR_<Direction::X>(J, nextX);

            auto JL_Y_1 = reconstruct_uL_<Direction::Y>(J, prevY);
            auto JR_Y_1 = reconstruct_uR_<Direction::Y>(J, prevY);
            auto JL_Y1  = reconstruct_uL_<Direction::Y>(J, nextY);
            auto JR_Y1  = reconstruct_uR_<Direction::Y>(J, nextY);

            auto JL_Z_1 = reconstruct_uL_<Direction::Z>(J, prevZ);
            auto JR_Z_1 = reconstruct_uR_<Direction::Z>(J, prevZ);
            auto JL_Z1  = reconstruct_uL_<Direction::Z>(J, nextZ);
            auto JR_Z1  = reconstruct_uR_<Direction::Z>(J, nextZ);

            auto LaplJL = d2(Direction::X, JL_X_1, JL, JL_X1) + d2(Direction::Y, JL_Y_1, JL, JL_Y1)
                          + d2(Direction::Z, JL_Z_1, JL, JL_Z1);
            auto LaplJR = d2(Direction::X, JR_X_1, JR, JR_X1) + d2(Direction::Y, JR_Y_1, JR, JR_Y1)
                          + d2(Direction::Z, JR_Z_1, JR, JR_Z1);

            return std::make_tuple(LaplJL, LaplJR);
        }
    }

    template<auto direction, typename VecField>
    auto center_reconstruct_(const VecField& U, MeshIndex<VecField::dimension> index,
                             auto projectionX, auto projectionY, auto projectionZ) const
    {
        auto const& Ux = U(Component::X);
        auto const& Uy = U(Component::Y);
        auto const& Uz = U(Component::Z);

        if (reconstruction_ == "constant")
        {
            auto UxL = GridLayout::project(Ux, previous_<direction>(Ux, index), projectionX);
            auto UyL = GridLayout::project(Uy, previous_<direction>(Uy, index), projectionY);
            auto UzL = GridLayout::project(Uz, previous_<direction>(Uz, index), projectionZ);

            auto UxR = GridLayout::project(Ux, index, projectionX);
            auto UyR = GridLayout::project(Uy, index, projectionY);
            auto UzR = GridLayout::project(Uz, index, projectionZ);

            return std::make_tuple(UxL, UyL, UzL, UxR, UyR, UzR);
        }
        else if (reconstruction_ == "linear")
        {
            auto UxL = center_recons_linear_L_<direction>(Ux, index, projectionX);
            auto UyL = center_recons_linear_L_<direction>(Uy, index, projectionY);
            auto UzL = center_recons_linear_L_<direction>(Uz, index, projectionZ);

            auto UxR = center_recons_linear_R_<direction>(Ux, index, projectionX);
            auto UyR = center_recons_linear_R_<direction>(Uy, index, projectionY);
            auto UzR = center_recons_linear_R_<direction>(Uz, index, projectionZ);

            return std::make_tuple(UxL, UyL, UzL, UxR, UyR, UzR);
        }
        else
            throw std::runtime_error("Error - Reconstruct - Unknown reconstruction method");
    }

    template<auto direction>
    auto center_recons_linear_L_(auto U, MeshIndex<dimension> index, auto projection) const
    {
        auto u_2 = GridLayout::project(U, previous_<direction>(U, previous_<direction>(U, index)),
                                       projection);
        auto u_1 = GridLayout::project(U, previous_<direction>(U, index), projection);
        auto u   = GridLayout::project(U, index, projection);

        return recons_linear_L_(u_2, u_1, u);
    }

    template<auto direction>
    auto center_recons_linear_R_(auto U, MeshIndex<dimension> index, auto projection) const
    {
        auto u_1 = GridLayout::project(U, previous_<direction>(U, index), projection);
        auto u   = GridLayout::project(U, index, projection);
        auto u1  = GridLayout::project(U, next_<direction>(U, index), projection);

        return recons_linear_R_(u_1, u, u1);
    }

    template<auto direction, typename Field>
    auto reconstruct_uL_(Field const& F, MeshIndex<Field::dimension> index) const
    {
        if (reconstruction_ == "constant")
            return constant_uL_<direction>(F, index);
        else if (reconstruction_ == "linear")
            return linear_uL_<direction>(F, index);
        else
            throw std::runtime_error("Error - Reconstruct - Unknown reconstruction method");
    }

    template<auto direction, typename Field>
    auto reconstruct_uR_(Field const& F, MeshIndex<Field::dimension> index) const
    {
        if (reconstruction_ == "constant")
            return constant_uR_<direction>(F, index);
        else if (reconstruction_ == "linear")
            return linear_uR_<direction>(F, index);
        else
            throw std::runtime_error("Error - Reconstruct - Unknown reconstruction method");
    }

    template<auto direction, typename Field>
    auto constant_uL_(Field const& F, MeshIndex<Field::dimension> index) const
    {
        return F(previous_<direction>(F, index));
    }

    template<auto direction, typename Field>
    auto constant_uR_(Field const& F, MeshIndex<Field::dimension> index) const
    {
        return F(index);
    }

    template<auto direction, typename Field>
    auto linear_uL_(Field const& F, MeshIndex<Field::dimension> index) const
    {
        auto ui_2 = F(previous_<direction>(F, previous_<direction>(F, index)));
        auto ui_1 = F(previous_<direction>(F, index));
        auto ui   = F(index);

        return recons_linear_L_(ui_2, ui_1, ui);
    }

    template<auto direction, typename Field>
    auto linear_uR_(Field const& F, MeshIndex<Field::dimension> index) const
    {
        auto ui_1 = F(previous_<direction>(F, index));
        auto ui   = F(index);
        auto ui1  = F(next_<direction>(F, index));

        return recons_linear_R_(ui_1, ui, ui1);
    }

    auto recons_linear_L_(auto ul, auto u, auto ur) const
    {
        auto Dil = (u - ul);
        auto Dir = (ur - u);
        auto Di  = slope_limiter_(Dil, Dir);

        return u + 0.5 * Di;
    }

    auto recons_linear_R_(auto ul, auto u, auto ur) const
    {
        auto Dil = (u - ul);
        auto Dir = (ur - u);
        auto Di  = slope_limiter_(Dil, Dir);

        return u - 0.5 * Di;
    }

    auto slope_limiter_(auto const& Dil, auto const& Dir) const
    {
        if (limiter_ == "minmod")
            return min_mod_limiter_(Dil, Dir);
        else if (limiter_ == "vanleer")
            return van_leer_limiter_(Dil, Dir);
        else
            throw std::runtime_error("Error - Slope Limiter - Unknown slope limiter method");
    }

    auto van_leer_limiter_(auto const& Dil, auto const& Dir) const
    {
        return Dil * Dir > 0.0 ? 2.0 * Dil * Dir / (Dil + Dir) : 0.0;
    }

    auto min_mod_limiter_(auto const& Dil, auto const& Dir) const
    {
        return Dil * Dir < 0.0 ? 0.0 : fabs(Dir) < fabs(Dil) ? Dir : Dil;
    }

    template<auto direction, typename Field>
    MeshIndex<Field::dimension> next_(Field const& F, MeshIndex<Field::dimension> index) const
    {
        auto fieldCentering = layout_->centering(F.physicalQuantity());

        if constexpr (dimension == 1)
        {
            return make_index(index[0] + 1);
        }
        else if constexpr (dimension == 2)
        {
            if constexpr (direction == Direction::X)
            {
                return make_index(index[0] + 1, index[1]);
            }
            else if constexpr (direction == Direction::Y)
            {
                return make_index(index[0], index[1] + 1);
            }
        }
        else if constexpr (dimension == 3)
        {
            if constexpr (direction == Direction::X)
            {
                return make_index(index[0] + 1, index[1], index[2]);
            }
            else if constexpr (direction == Direction::Y)
            {
                return make_index(index[0], index[1] + 1, index[2]);
            }
            else if constexpr (direction == Direction::Z)
            {
                return make_index(index[0], index[1], index[2] + 1);
            }
        }
    }

    template<auto direction, typename Field>
    MeshIndex<Field::dimension> previous_(Field const& F, MeshIndex<Field::dimension> index) const
    {
        auto fieldCentering = layout_->centering(F.physicalQuantity());

        if constexpr (dimension == 1)
        {
            return make_index(index[0] - 1);
        }
        else if constexpr (dimension == 2)
        {
            if constexpr (direction == Direction::X)
            {
                return make_index(index[0] - 1, index[1]);
            }
            else if constexpr (direction == Direction::Y)
            {
                return make_index(index[0], index[1] - 1);
            }
        }
        else if constexpr (dimension == 3)
        {
            if constexpr (direction == Direction::X)
            {
                return make_index(index[0] - 1, index[1], index[2]);
            }
            else if constexpr (direction == Direction::Y)
            {
                return make_index(index[0], index[1] - 1, index[2]);
            }
            else if constexpr (direction == Direction::Z)
            {
                return make_index(index[0], index[1], index[2] - 1);
            }
        }
    }
};

} // namespace PHARE::core

#endif
