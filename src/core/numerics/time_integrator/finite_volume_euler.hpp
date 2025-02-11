#ifndef PHARE_CORE_NUMERICS_EULER_HPP
#define PHARE_CORE_NUMERICS_EULER_HPP

#include "initializer/data_provider.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler_per_field.hpp"

namespace PHARE::core
{
template<typename GridLayout>
class FiniteVolumeEuler : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    FiniteVolumeEuler(PHARE::initializer::PHAREDict const& dict)
        : gamma_{dict["heat_capacity_ratio"].template to<double>()}
    {
    }

    template<typename State, typename... Fluxes>
    void operator()(State const& state, State& statenew, double const dt,
                    Fluxes const&... fluxes) const
    {
        auto const fve = FiniteVolumeEulerPerField_ref{*layout_, dt};

        auto& rhoVxnew = statenew.rhoV(Component::X);
        auto& rhoVynew = statenew.rhoV(Component::Y);
        auto& rhoVznew = statenew.rhoV(Component::Z);

        auto const& rhoVx = state.rhoV(Component::X);
        auto const& rhoVy = state.rhoV(Component::Y);
        auto const& rhoVz = state.rhoV(Component::Z);

        auto const& f = std::forward_as_tuple(fluxes...);

        auto const& rho_fx  = std::get<0>(f);
        auto const& rhoV_fx = std::get<1>(f);
        auto const& Etot_fx = std::get<3>(f);

        auto const& rhoVx_fx = rhoV_fx(Component::X);
        auto const& rhoVy_fx = rhoV_fx(Component::Y);
        auto const& rhoVz_fx = rhoV_fx(Component::Z);

        if constexpr (dimension == 1)
        {
            fve(state.rho, statenew.rho, rho_fx);
            fve(rhoVx, rhoVxnew, rhoVx_fx);
            fve(rhoVy, rhoVynew, rhoVy_fx);
            fve(rhoVz, rhoVznew, rhoVz_fx);
            fve(state.Etot, statenew.Etot, Etot_fx);
        }

        if constexpr (dimension >= 2)
        {
            auto const& rho_fy  = std::get<4>(f);
            auto const& rhoV_fy = std::get<5>(f);
            auto const& Etot_fy = std::get<7>(f);

            auto const& rhoVx_fy = rhoV_fy(Component::X);
            auto const& rhoVy_fy = rhoV_fy(Component::Y);
            auto const& rhoVz_fy = rhoV_fy(Component::Z);

            if constexpr (dimension == 2)
            {
                fve(state.rho, statenew.rho, rho_fx, rho_fy);
                fve(rhoVx, rhoVxnew, rhoVx_fx, rhoVx_fy);
                fve(rhoVy, rhoVynew, rhoVy_fx, rhoVy_fy);
                fve(rhoVz, rhoVznew, rhoVz_fx, rhoVz_fy);
                fve(state.Etot, statenew.Etot, Etot_fx, Etot_fy);
            }
            if constexpr (dimension == 3)
            {
                auto const& rho_fz  = std::get<8>(f);
                auto const& rhoV_fz = std::get<9>(f);
                auto const& Etot_fz = std::get<11>(f);

                auto const& rhoVx_fz = rhoV_fz(Component::X);
                auto const& rhoVy_fz = rhoV_fz(Component::Y);
                auto const& rhoVz_fz = rhoV_fz(Component::Z);

                fve(state.rho, statenew.rho, rho_fx, rho_fy, rho_fz);
                fve(rhoVx, rhoVxnew, rhoVx_fx, rhoVx_fy, rhoVx_fz);
                fve(rhoVy, rhoVynew, rhoVy_fx, rhoVy_fy, rhoVy_fz);
                fve(rhoVz, rhoVznew, rhoVz_fx, rhoVz_fy, rhoVz_fz);
                fve(state.Etot, statenew.Etot, Etot_fx, Etot_fy, Etot_fz);
            }
        }
    }

private:
    double const gamma_;
};

} // namespace PHARE::core

#endif
