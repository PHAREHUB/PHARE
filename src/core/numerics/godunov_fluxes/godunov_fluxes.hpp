#ifndef PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP
#define PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/point/point.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/reconstructions/reconstructor.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/types.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <tuple>

namespace PHARE::core
{
template<size_t dim>
constexpr auto getDirections()
{
    if constexpr (dim == 1)
    {
        return std::make_tuple(Direction::X);
    }
    else if constexpr (dim == 2)
    {
        return std::make_tuple(Direction::X, Direction::Y);
    }
    else if constexpr (dim == 3)
    {
        return std::make_tuple(Direction::X, Direction::Y, Direction::Z);
    }
}

template<auto direction, size_t dim, bool HyperResistivity>
auto getGrow(int const nghosts)
{
    Point<std::uint32_t, dim> p{};

    auto dir = static_cast<size_t>(direction);

    for (size_t i = 0; i < dim; ++i)
    {
        if (i != dir)
            p[i] = nghosts;
    }

    // add one extra layer in the direction of the flux laplacian computation. Maybe some later
    // optimisation would let us just compute for uct and have the extra layer only reconstructed
    // for j
    if constexpr (HyperResistivity)
    {
        p[dir] += 1;
    }

    return p;
}

template<typename GridLayout, typename MHDModel, template<typename> typename Reconstruction,
         typename RiemannSolver, typename Equations>
class Godunov : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

    using Reconstruction_t = Reconstruction<GridLayout>;
    using Reconstructor_t  = Reconstructor<Reconstruction_t>;
    using RiemannSolver_t  = RiemannSolver;

public:
    template<typename T>
    using Rec = Reconstruction<T>;

    constexpr static auto Hall             = Equations::hall;
    constexpr static auto Resistivity      = Equations::resistivity;
    constexpr static auto HyperResistivity = Equations::hyperResistivity;

    Godunov(PHARE::initializer::PHAREDict const& dict)
        : gamma_{dict["heat_capacity_ratio"].template to<double>()}
        , eta_{dict["resistivity"].template to<double>()}
        , nu_{dict["hyper_resistivity"].template to<double>()}
        , hyper_mode_{cppdict::get_value(dict, "hyper_mode", std::string{"constant"}) == "constant"
                          ? HyperMode::constant
                          : HyperMode::spatial}
        , equations_{gamma_, eta_, nu_}
        , riemann_{gamma_}
    {
    }

    template<typename State, typename Fluxes>
    void operator()(auto& ct, State& state, Fluxes& fluxes)
    {
        if (!this->hasLayout())
            throw std::runtime_error("Error - GodunovFluxes - GridLayout not set");

        constexpr auto directions = getDirections<dimension>();

        constexpr auto num_directions = std::tuple_size_v<std::decay_t<decltype(directions)>>;

        for_N<num_directions>([&](auto i) {
            constexpr Direction direction = std::get<i>(directions);

            layout_->evalOnBiggerBox(
                fluxes.template expose_centering<direction>(),
                getGrow<direction, dimension, HyperResistivity>(Reconstruction_t::nghosts),
                [&](auto&... indices) {
                    if constexpr (Hall || Resistivity || HyperResistivity)
                    {
                        auto&& [uL, uR]
                            = Reconstructor_t::template reconstruct<direction>(state, {indices...});

                        auto const& [jL, jR]
                            = Reconstructor_t::template center_reconstruct<
                                direction, GridLayout::edgeXToCellCenter,
                                GridLayout::edgeYToCellCenter, GridLayout::edgeZToCellCenter>(
                                state.J, {indices...});

                        auto&& u      = std::forward_as_tuple(uL, uR);
                        auto const& j = std::forward_as_tuple(jL, jR);

                        // if constexpr (HyperResistivity)
                        // {
                        //     auto const& [laplJL, laplJR]
                        //         = Reconstructor_t::template reconstructed_laplacian<direction>(
                        //             layout_->inverseMeshSize(), state.J, {indices...});
                        //
                        //     auto const& LaplJ = std::forward_as_tuple(laplJL, laplJR);
                        //
                        //     auto const& [fL, fR] = for_N<2, for_N_R_mode::make_tuple>([&](auto i)
                        //     {
                        //         return equations_.template compute<direction>(
                        //             std::get<i>(u), std::get<i>(j), std::get<i>(LaplJ));
                        //     });
                        //
                        //     fluxes.template get_dir<direction>({indices...})
                        //         = riemann_.template solve<direction>(uL, uR, fL, fR, jL, jR);
                        //
                        //     ct.template save<direction>(riemann_.vt, riemann_.jt,
                        //                                 riemann_.rhot, riemann_.uct_coefs,
                        //                                 {indices...});
                        // }
                        // else
                        // {
                        auto const& [fL, fR] = for_N<2, for_N_R_mode::make_tuple>([&](auto i) {
                            return equations_.template compute<direction>(std::get<i>(u),
                                                                          std::get<i>(j));
                        });

                        // if constexpr (Hall)
                        // {
                        fluxes.template get_dir<direction>({indices...})
                            = riemann_.template solve<direction>(uL, uR, fL, fR, jL, jR);

                        ct.template save<direction>(riemann_.vt, riemann_.jt, riemann_.rhot,
                                                    riemann_.uct_coefs, {indices...});
                        // }
                        // else // Resistivity only
                        // {
                        //     fluxes.template get_dir<direction>({indices...})
                        //         = riemann_.template solve<direction>(uL, uR, fL, fR);
                        //
                        //     ct.template save<direction>(riemann_.vt,
                        //                                 riemann_.uct_coefs, {indices...});
                        // }
                        // }
                    }
                    else // Ideal
                    {
                        auto&& [uL, uR]
                            = Reconstructor_t::template reconstruct<direction>(state, {indices...});

                        auto&& u = std::forward_as_tuple(uL, uR);

                        auto const& [fL, fR] = for_N<2, for_N_R_mode::make_tuple>([&](auto i) {
                            return equations_.template compute<direction>(std::get<i>(u));
                        });

                        fluxes.template get_dir<direction>({indices...})
                            = riemann_.template solve<direction>(uL, uR, fL, fR);

                        ct.template save<direction>(riemann_.vt, riemann_.uct_coefs, {indices...});
                    }
                });

            // Note: resistive contributions to F_B and F_Etot are now handled via CT:
            // - B is updated by Faraday using E (which includes ηJ from CT)
            // - Etot gets the resistive Poynting flux via E×B where E = E_ideal + ηJ
            // The old resistive_contributions code was dead (F_B not used) and would
            // double-count η(J×B) in F_Etot when combined with Poynting correction.
        });
    }

    void registerResources(MHDModel& model) {}

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const {}

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(); }

    template<typename CT, typename State, typename Fluxes>
    void apply_poynting_correction(CT const& ct, State const& state, Fluxes& fluxes)
    {
        // Apply Poynting flux correction to energy: ∂E/∂t -= ∇·(E×B)
        // This must be called AFTER CT has computed both E and edge-B fields
        // Uses edge-centered B from CT (guaranteed temporally consistent with E)
        //
        // The E field from CT includes resistive contributions (E = -V×B + ηJ + ...),
        // so E×B automatically captures the resistive Poynting flux η(J×B).
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - GodunovFluxes::apply_poynting_correction - GridLayout not set");

        constexpr auto directions     = getDirections<dimension>();
        constexpr auto num_directions = std::tuple_size_v<std::decay_t<decltype(directions)>>;

        for_N<num_directions>([&](auto i) {
            constexpr Direction direction = std::get<i>(directions);

            layout_->evalOnBox(
                fluxes.template expose_centering<direction>(), [&](auto&... indices) {
                    auto& F_Etot = fluxes.template get_dir<direction>({indices...}).Etot();
                    poynting_energy_flux_<direction>(ct, state.E, MeshIndex<dimension>{indices...},
                                                     F_Etot);
                });
        });
    }

private:
    template<auto direction>
    void poynting_energy_flux_(auto const& ct, auto const& E, MeshIndex<dimension> const& index,
                               auto& F_Etot) const
    {
        // Compute magnetic energy flux via Poynting vector: S·n̂ = (E × B)·n̂
        // E components live on edges (from CT)
        // B components are edge-centered (from CT, temporally consistent with E)
        //
        // gPluto pattern: average E×B products perpendicular to the edge direction
        // to bring the product from edge-centered to face-centered

        auto const& Ex = E(Component::X);
        auto const& Ey = E(Component::Y);
        auto const& Ez = E(Component::Z);

        if constexpr (direction == Direction::X && dimension >= 2)
        {
            // X-flux face: Sx = Ey*Bz - Ez*By
            // Ez lives at z-edge, By_at_Ez lives at z-edge -> average in y-direction
            // Ey lives at y-edge (= x-face in 2D), Bz_at_Ey at same location
            //   In 3D: average in z-direction; in 2D: no averaging needed (same face)

            auto const& By_at_Ez = ct.getBy_at_Ez(); // gPluto: Bx2ez
            auto const& Bz_at_Ey = ct.getBz_at_Ey(); // gPluto: Bx3ey (2D+)

            double EzBy = 0.5
                          * (Ez(index) * By_at_Ez(index)
                             + Ez(layout_->template next<Direction::Y>(index))
                                   * By_at_Ez(layout_->template next<Direction::Y>(index)));

            double EyBz;
            if constexpr (dimension == 2)
            {
                // In 2D, Ey lives at y-edges = x-faces: same location as x-flux face.
                // No averaging needed.
                EyBz = Ey(index) * Bz_at_Ey(index);
            }
            else // dimension == 3
            {
                // gPluto: EyBz = 0.5*(Ey(i,j,k)*Bz_y(i,j,k) + Ey(i,j,k+1)*Bz_y(i,j,k+1))
                EyBz = 0.5
                       * (Ey(index) * Bz_at_Ey(index)
                          + Ey(layout_->template next<Direction::Z>(index))
                                * Bz_at_Ey(layout_->template next<Direction::Z>(index)));
            }

            F_Etot += EyBz - EzBy;
        }
        else if constexpr (direction == Direction::Y && dimension >= 2)
        {
            // Y-flux face: Sy = Ez*Bx - Ex*Bz
            // Ez lives at z-edge, Bx_at_Ez lives at z-edge -> average in x-direction
            // Ex lives at x-edge (= y-face in 2D), Bz_at_Ex at same location
            //   In 3D: average in z-direction; in 2D: no averaging needed

            auto const& Bx_at_Ez = ct.getBx_at_Ez(); // gPluto: Bx1ez
            auto const& Bz_at_Ex = ct.getBz_at_Ex(); // gPluto: Bx3ex (2D+)

            double EzBx = 0.5
                          * (Ez(index) * Bx_at_Ez(index)
                             + Ez(layout_->template next<Direction::X>(index))
                                   * Bx_at_Ez(layout_->template next<Direction::X>(index)));

            double ExBz;
            if constexpr (dimension == 2)
            {
                // In 2D, Ex lives at x-edges = y-faces: same location as y-flux face.
                // No averaging needed.
                ExBz = Ex(index) * Bz_at_Ex(index);
            }
            else // dimension == 3
            {
                // gPluto: ExBz = 0.5*(Ex(i,j,k)*Bz_x(i,j,k) + Ex(i,j,k+1)*Bz_x(i,j,k+1))
                ExBz = 0.5
                       * (Ex(index) * Bz_at_Ex(index)
                          + Ex(layout_->template next<Direction::Z>(index))
                                * Bz_at_Ex(layout_->template next<Direction::Z>(index)));
            }

            F_Etot += EzBx - ExBz;
        }
        else if constexpr (direction == Direction::Z && dimension == 3)
        {
            // Z-flux face: Sz = Ex*By - Ey*Bx
            // Ex lives at x-edge, By_at_Ex lives at x-edge -> average in y-direction
            // Ey lives at y-edge, Bx_at_Ey lives at y-edge -> average in x-direction

            auto const& By_at_Ex = ct.getBy_at_Ex(); // gPluto: Bx2ex
            auto const& Bx_at_Ey = ct.getBx_at_Ey(); // gPluto: Bx1ey

            // gPluto: ExBy = 0.5*(Ex(i,j,k)*By_x(i,j,k) + Ex(i,j+1,k)*By_x(i,j+1,k))
            double ExBy = 0.5
                          * (Ex(index) * By_at_Ex(index)
                             + Ex(layout_->template next<Direction::Y>(index))
                                   * By_at_Ex(layout_->template next<Direction::Y>(index)));

            // gPluto: EyBx = 0.5*(Ey(i,j,k)*Bx_y(i,j,k) + Ey(i+1,j,k)*Bx_y(i+1,j,k))
            double EyBx = 0.5
                          * (Ey(index) * Bx_at_Ey(index)
                             + Ey(layout_->template next<Direction::X>(index))
                                   * Bx_at_Ey(layout_->template next<Direction::X>(index)));

            F_Etot += ExBy - EyBx;
        }
        // direction == X && dimension == 1: No Poynting correction (no transverse directions)
    }


    double const gamma_;
    double const eta_;
    double const nu_;
    HyperMode const hyper_mode_;

    Equations equations_;
    RiemannSolver_t riemann_;
};

} // namespace PHARE::core

#endif
