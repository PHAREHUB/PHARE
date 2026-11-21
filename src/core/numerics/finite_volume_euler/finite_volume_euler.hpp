#ifndef PHARE_CORE_NUMERICS_EULER_HPP
#define PHARE_CORE_NUMERICS_EULER_HPP


#include "core/utilities/types.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler_per_field.hpp"

#include <tuple>
#include <utility>

namespace PHARE::core
{
template<typename GridLayout>
class FiniteVolumeEuler
{
    constexpr static auto dimension = GridLayout::dimension;


public:
    FiniteVolumeEuler(GridLayout const& layout)
        : layout_{layout}

    {
    }

    template<typename State, typename Fluxes>
    void operator()(State const& state, State& statenew, Fluxes const& fluxes,
                    double const dt) const
    {
        auto const fve = FiniteVolumeEulerPerField{layout_, dt};
        update_(state, statenew, [&](auto const& U, auto& Unew, auto const&... flux) {
            fve(U, Unew, flux(fluxes)...);
        });
    }

    // the update with the fluxes A + w B, their sum formed face by face instead of stored
    template<typename State, typename Fluxes>
    void operator()(State const& state, State& statenew, Fluxes const& A, double const w,
                    Fluxes const& B, double const dt) const
    {
        auto const fve = FiniteVolumeEulerPerField{layout_, dt};
        update_(state, statenew, [&](auto const& U, auto& Unew, auto const&... flux) {
            fve.weighted(U, Unew, w,
                         std::pair<decltype(flux(A)), decltype(flux(B))>{flux(A), flux(B)}...);
        });
    }

private:
    // fn(U, Unew, flux_x[, flux_y[, flux_z]]) for each conserved field U, flux_d(fluxes) being
    // its flux along direction d
    template<typename State>
    void update_(State const& state, State& statenew, auto&& fn) const
    {
        auto const apply = [&](auto const& U, auto& Unew, auto const& select) {
            [&]<std::size_t... d>(std::index_sequence<d...>) {
                fn(U, Unew, [&](auto const& F) -> auto const& {
                    return select(
                        std::get<d>(std::forward_as_tuple(F.rho_fx, F.rho_fy, F.rho_fz)),
                        std::get<d>(std::forward_as_tuple(F.rhoV_fx, F.rhoV_fy, F.rhoV_fz)),
                        std::get<d>(std::forward_as_tuple(F.Etot_fx, F.Etot_fy, F.Etot_fz)));
                }...);
            }(std::make_index_sequence<dimension>{});
        };

        apply(state.rho, statenew.rho,
              [](auto const& rho, auto const&, auto const&) -> auto const& { return rho; });
        for_N<3>([&](auto c) {
            constexpr auto component = static_cast<Component>(c());
            apply(state.rhoV[component], statenew.rhoV[component],
                  [component](auto const&, auto const& rhoV, auto const&) -> auto const& {
                      return rhoV[component];
                  });
        });
        apply(state.Etot, statenew.Etot,
              [](auto const&, auto const&, auto const& Etot) -> auto const& { return Etot; });
    }

    GridLayout layout_;
};


} // namespace PHARE::core

#endif
