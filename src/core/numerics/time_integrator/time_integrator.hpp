#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_HPP


#include "core/data/grid/gridlayout_utils.hpp"
#include "core/numerics/constrained_transport/constrained_transport.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler.hpp"
#include "core/numerics/godunov_fluxes/godunov_fluxes.hpp"
#include "initializer/data_provider.hpp"
#include <string>

namespace PHARE::core
{
template<typename GridLayout>
class TimeIntegrator : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename Field, typename VecField, typename... Fluxes>
    void euler(Field& rho, VecField& rhoV, VecField& B, Field& Etot, VecField& E, double const dt,
               Fluxes&... fluxes) const
    {
        euler_step_(rho, rhoV, B, Etot, rho, rhoV, B, Etot, E, dt, fluxes...);
    }

    template<typename Field, typename VecField, typename... Fluxes>
    void tvdrk2(Field& rho, VecField& rhoV, VecField& B, Field& Etot, Field& rho1, VecField& rhoV1,
                VecField& B1, Field& Etot1, VecField& E, double const dt, Fluxes&... fluxes) const
    {
        euler_step_(rho, rhoV, B, Etot, rho1, rhoV1, B1, Etot1, E, dt,
                    fluxes...); // step 1: U1 = Euler(Un)

        euler_step_(rho1, rhoV1, B1, Etot1, rho1, rhoV1, B1, Etot1, E, dt,
                    fluxes...); // U1 <- Euler(U1)

        auto& rhoVx = rhoV(Component::X);
        auto& rhoVy = rhoV(Component::Y);
        auto& rhoVz = rhoV(Component::Z);

        auto& Bx = B(Component::X);
        auto& By = B(Component::Y);
        auto& Bz = B(Component::Z);

        auto& rhoVx1 = rhoV1(Component::X);
        auto& rhoVy1 = rhoV1(Component::Y);
        auto& rhoVz1 = rhoV1(Component::Z);

        auto& Bx1 = B1(Component::X);
        auto& By1 = B1(Component::Y);
        auto& Bz1 = B1(Component::Z);

        auto U  = std::forward_as_tuple(rho, rhoVx, rhoVy, rhoVz, Bx, By, Bz, Etot);
        auto U1 = std::forward_as_tuple(rho1, rhoVx1, rhoVy1, rhoVz1, Bx1, By1, Bz1, Etot1);

        auto constexpr N_elements = std::tuple_size_v<std::decay_t<decltype(U)>>;

        auto step_2 // step 2: Un1 = 0.5 Un + 0.5 Euler(U1)
            = [&](Field& field, Field& field1, MeshIndex<dimension> index) {
                  field(index) = 0.5 * field(index) + 0.5 * field1(index);
              };

        for_N<N_elements>([&](auto i) {
            layout_->evalOnBox(std::get<i>(U), [&](auto... indices) mutable {
                step_2(std::get<i>(U), std::get<i>(U1), {indices...});
            });
        });
    }

    template<typename Field, typename VecField, typename... Fluxes>
    void tvdrk3(Field& rho, VecField& rhoV, VecField& B, Field& Etot, Field& rho1, VecField& rhoV1,
                VecField& B1, Field& Etot1, Field& rho2, VecField& rhoV2, VecField& B2,
                Field& Etot2, VecField& E, double const dt, Fluxes&... fluxes) const
    {
        euler_step_(rho, rhoV, B, Etot, rho1, rhoV1, B1, Etot1, E, dt,
                    fluxes...); // step 1: U1 = Euler(Un)

        euler_step_(rho1, rhoV1, B1, Etot1, rho1, rhoV1, B1, Etot1, E, dt,
                    fluxes...); // U1 <- Euler(U1)

        auto& rhoVx = rhoV(Component::X);
        auto& rhoVy = rhoV(Component::Y);
        auto& rhoVz = rhoV(Component::Z);

        auto& Bx = B(Component::X);
        auto& By = B(Component::Y);
        auto& Bz = B(Component::Z);

        auto& rhoVx1 = rhoV1(Component::X);
        auto& rhoVy1 = rhoV1(Component::Y);
        auto& rhoVz1 = rhoV1(Component::Z);

        auto& Bx1 = B1(Component::X);
        auto& By1 = B1(Component::Y);
        auto& Bz1 = B1(Component::Z);

        auto& rhoVx2 = rhoV2(Component::X);
        auto& rhoVy2 = rhoV2(Component::Y);
        auto& rhoVz2 = rhoV2(Component::Z);

        auto& Bx2 = B2(Component::X);
        auto& By2 = B2(Component::Y);
        auto& Bz2 = B2(Component::Z);

        auto U  = std::forward_as_tuple(rho, rhoVx, rhoVy, rhoVz, Bx, By, Bz, Etot);
        auto U1 = std::forward_as_tuple(rho1, rhoVx1, rhoVy1, rhoVz1, Bx1, By1, Bz1, Etot1);
        auto U2 = std::forward_as_tuple(rho2, rhoVx2, rhoVy2, rhoVz2, Bx2, By2, Bz2, Etot2);

        auto constexpr N_elements = std::tuple_size_v<std::decay_t<decltype(U)>>;

        auto step_2 // step 2: U2 = 0.75 Un + 0.25 Euler(U1)
            = [&](auto& field, auto& field1, auto& field2, MeshIndex<dimension> index) {
                  field2(index) = 0.75 * field(index) + 0.25 * field1(index);
              };

        for_N<N_elements>([&](auto i) {
            layout_->evalOnBox(std::get<i>(U2), [&](auto... indices) mutable {
                step_2(std::get<i>(U), std::get<i>(U1), std::get<i>(U2), {indices...});
            });
        });

        euler_step_(rho2, rhoV2, B2, Etot2, rho2, rhoV2, B2, Etot2, E, dt,
                    fluxes...); // U2 <- Euler(U2)

        auto step_3 // step 3: Un1 = 1/3 Un + 2/3 Euler(U2)
            = [&](Field& field, Field& field2, MeshIndex<dimension> index) {
                  field(index) = 1.0 / 3.0 * field(index) + 2.0 / 3.0 * field2(index);
              };

        for_N<N_elements>([&](auto i) {
            layout_->evalOnBox(std::get<i>(U), [&](auto... indices) mutable {
                step_3(std::get<i>(U), std::get<i>(U2), {indices...});
            });
        });
    }

private:
    template<typename Field, typename VecField, typename... Fluxes>
    void euler_step_(Field const& rho, VecField const& rhoV, VecField const& B, Field const& Etot,
                     Field& rhonew, VecField& rhoVnew, VecField& Bnew, Field& Etotnew, VecField& E,
                     double const dt, Fluxes const&... fluxes) const
    {
        auto const fve     = FiniteVolumeEuler_ref{*layout_, dt};
        auto const ct      = ConstrainedTransport_ref{*layout_};
        auto const faraday = Faraday_ref{*layout_, dt};

        auto&& flux_tuple = std::forward_as_tuple(fluxes...);

        auto& rho_x  = std::get<0>(flux_tuple);
        auto& rhoV_x = std::get<1>(flux_tuple);
        auto& B_x    = std::get<2>(flux_tuple);
        auto& Etot_x = std::get<3>(flux_tuple);

        auto& rhoVxnew = rhoVnew(Component::X);
        auto& rhoVynew = rhoVnew(Component::Y);
        auto& rhoVznew = rhoVnew(Component::Z);

        auto const& rhoVx = rhoV(Component::X);
        auto const& rhoVy = rhoV(Component::Y);
        auto const& rhoVz = rhoV(Component::Z);

        auto const& rhoVx_x = rhoV_x(Component::X);
        auto const& rhoVy_x = rhoV_x(Component::Y);
        auto const& rhoVz_x = rhoV_x(Component::Z);


        if constexpr (dimension == 1)
        {
            fve(rho, rhonew, rho_x);
            fve(rhoVx, rhoVxnew, rhoVx_x);
            fve(rhoVy, rhoVynew, rhoVy_x);
            fve(rhoVz, rhoVznew, rhoVz_x);
            fve(Etot, Etotnew, Etot_x);

            ct(E, B_x);
            faraday.onBox(B, E, Bnew);
        }

        if constexpr (dimension >= 2)
        {
            auto& rho_y  = std::get<4>(flux_tuple);
            auto& rhoV_y = std::get<5>(flux_tuple);
            auto& B_y    = std::get<6>(flux_tuple);
            auto& Etot_y = std::get<7>(flux_tuple);

            auto const& rhoVx_y = rhoV_y(Component::X);
            auto const& rhoVy_y = rhoV_y(Component::Y);
            auto const& rhoVz_y = rhoV_y(Component::Z);

            if constexpr (dimension == 2)
            {
                fve(rho, rhonew, rho_x, rho_y);
                fve(rhoVx, rhoVxnew, rhoVx_x, rhoVx_y);
                fve(rhoVy, rhoVynew, rhoVy_x, rhoVy_y);
                fve(rhoVz, rhoVznew, rhoVz_x, rhoVz_y);
                fve(Etot, Etotnew, Etot_x, Etot_y);

                ct(E, B_x, B_y);
                faraday.onBox(B, E, Bnew);
            }
            if constexpr (dimension == 3)
            {
                auto& rho_z  = std::get<8>(flux_tuple);
                auto& rhoV_z = std::get<9>(flux_tuple);
                auto& B_z    = std::get<10>(flux_tuple);
                auto& Etot_z = std::get<11>(flux_tuple);

                auto const& rhoVx_z = rhoV_z(Component::X);
                auto const& rhoVy_z = rhoV_z(Component::Y);
                auto const& rhoVz_z = rhoV_z(Component::Z);

                fve(rho, rhonew, rho_x, rho_y, rho_z);
                fve(rhoVx, rhoVxnew, rhoVx_x, rhoVx_y, rhoVx_z);
                fve(rhoVy, rhoVynew, rhoVy_x, rhoVy_y, rhoVy_z);
                fve(rhoVz, rhoVznew, rhoVz_x, rhoVz_y, rhoVz_z);
                fve(Etot, Etotnew, Etot_x, Etot_y, Etot_z);

                ct(E, B_x, B_y, B_z);
                faraday.onBox(B, E, Bnew);
            }
        }
    }
};

} // namespace PHARE::core


#endif
