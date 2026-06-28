#ifndef PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP
#define PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/reconstructions/reconstructor.hpp"

#include "initializer/data_provider.hpp"

#include <cmath>
#include <tuple>
#include <cstddef>
#include <cstdint>

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
        if (i != dir)
            p[i] = nghosts;

    // add one extra layer in the direction of the flux laplacian computation. Maybe some later
    // optimisation would let us just compute for uct and have the extra layer only reconstructed
    // for j
    if constexpr (HyperResistivity)
        p[dir] += 1;

    return p;
}

struct GodunovInfo : public OhmInfo
{
    double const gamma;

    GodunovInfo static FROM(initializer::PHAREDict const& dict)
    {
        return {{OhmInfo::FROM(dict)}, dict["heat_capacity_ratio"].template to<double>()};
    }
};


template<typename GridLayout, template<typename> typename Reconstruction, typename RiemannSolver,
         typename Equations>
class Godunov : public GodunovInfo
{
    using Super                     = GodunovInfo;
    using Reconstruction_t          = Reconstruction<GridLayout>;
    using Reconstructor_t           = Reconstructor<Reconstruction_t>;
    using RiemannSolver_t           = RiemannSolver;
    constexpr static auto dimension = GridLayout::dimension;

public:
    using Info_t      = Super;
    using Equations_t = Equations;

    template<typename T>
    using Rec = Reconstruction<T>;

    constexpr static auto Hall             = Equations::hall;
    constexpr static auto Resistivity      = Equations::resistivity;
    constexpr static auto HyperResistivity = Equations::hyperResistivity;

    explicit Godunov(GodunovInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
        , equations_{gamma, eta, nu}
        , riemann_{gamma}
    {
    }

    template<typename State, typename Fluxes>
    void operator()(auto& ct_state, State& state, Fluxes& fluxes)
    {
        constexpr auto directions = getDirections<dimension>();

        constexpr auto num_directions = std::tuple_size_v<std::decay_t<decltype(directions)>>;

        for_N<num_directions>([&](auto i) {
            constexpr Direction direction = std::get<i>(directions);

            layout_.evalOnBiggerBox(
                fluxes.template expose_centering<direction>(),
                getGrow<direction, dimension, HyperResistivity>(Reconstruction_t::nghosts),
                [&](auto&... indices) {
                    if constexpr (Hall || Resistivity || HyperResistivity)
                    {
                        auto&& [uL, uR]
                            = Reconstructor_t::template reconstruct<direction>(state, {indices...});

                        auto const& [jL, jR] = Reconstructor_t::template center_reconstruct<
                            direction, GridLayout::implT::edgeXToCellCenter,
                            GridLayout::implT::edgeYToCellCenter,
                            GridLayout::implT::edgeZToCellCenter>(state.J, {indices...});

                        auto&& u      = std::forward_as_tuple(uL, uR);
                        auto const& j = std::forward_as_tuple(jL, jR);


                        auto const& [fL, fR] = for_N<2, for_N_R_mode::make_tuple>([&](auto i) {
                            return equations_.template compute<direction>(std::get<i>(u),
                                                                          std::get<i>(j));
                        });

                        fluxes.template get_dir<direction>({indices...})
                            = riemann_.template solve<direction>(uL, uR, fL, fR, jL, jR);

                        ct_state.template save<direction>(riemann_.vt, riemann_.jt, riemann_.rhot,
                                                          riemann_.uct_coefs, {indices...});
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

                        ct_state.template save<direction>(riemann_.vt, riemann_.uct_coefs,
                                                          {indices...});
                    }
                });

            // Note: resistive contributions to F_B and F_Etot are now handled via CT:
            // - B is updated by Faraday using E (which includes ηJ from CT)
            // - Etot gets the resistive Poynting flux via E×B where E = E_ideal + ηJ
            // The old resistive_contributions code was dead (F_B is never consumed by the
            // finite-volume update) and would double-count η(J×B) in F_Etot once combined
            // with the Poynting correction.
        });
    }

    template<typename CTState, typename State, typename Fluxes>
    void apply_poynting_correction(CTState const& ct_state, State const& state, Fluxes& fluxes)
    {
        // Apply Poynting flux correction to energy: ∂E/∂t -= ∇·(E×B)
        // Must be called AFTER CT has computed both E and the edge-B fields, so the
        // edge-centered B from CT is temporally consistent with E. Because E from CT
        // includes resistive (ηJ) and hyper-resistive terms, E×B automatically captures
        // the resistive Poynting flux.
        constexpr auto directions     = getDirections<dimension>();
        constexpr auto num_directions = std::tuple_size_v<std::decay_t<decltype(directions)>>;

        for_N<num_directions>([&](auto i) {
            constexpr Direction direction = std::get<i>(directions);

            layout_.evalOnBox(fluxes.template expose_centering<direction>(), [&](auto&... indices) {
                auto& F_Etot = fluxes.template get_dir<direction>({indices...}).Etot();
                poynting_energy_flux_<direction>(ct_state, state.E,
                                                 MeshIndex<dimension>{indices...}, F_Etot);
            });
        });
    }

private:
    template<auto direction>
    void poynting_energy_flux_(auto const& ct_state, auto const& E,
                               MeshIndex<dimension> const& index, auto& F_Etot) const
    {
        // Magnetic energy flux through the face: S.n = (E x B).n. Since S_i = eps_ijk E_j B_k,
        // a face never needs its own E component, which leaves exactly two products per face.
        // E and the CT edge-B buffers are point values living on edges; each product is formed
        // on the edge it lives on, then projected onto the flux face across the single axis
        // separating the two locations (primal -> dual in every case).
        //
        // directionalInterp returns the identity stencil for an axis the simulation does not
        // have, which is what collapses the lower dimensional cases: in 2D Ey already sits on
        // the x-face, and in 1D both products already sit on the x-face.
        using InterpDir = typename GridLayout::implT::InterpDir;

        // Same as GridLayout::project, but for an expression rather than a stored field: the
        // integrand is a product of two fields, and since both are point values the product
        // must be formed on the edge before being projected onto the face. Projecting the
        // factors separately and multiplying on the face would be a different (and, beyond
        // order 2, wrong) discretisation.
        auto project_on_face = [&](auto const& stencil, auto&& integrand) {
            std::decay_t<decltype(integrand(index))> result = 0.;

            for (auto const& wp : stencil)
                result += wp.coef * integrand(index + wp.indexes);

            return result;
        };

        auto const& Ex = E(Component::X);
        auto const& Ey = E(Component::Y);
        auto const& Ez = E(Component::Z);

        if constexpr (direction == Direction::X)
        {
            constexpr auto zEdgeToXFace
                = GridLayout::implT::template directionalInterp<dirY, InterpDir::PrimalToDual>();
            constexpr auto yEdgeToXFace
                = GridLayout::implT::template directionalInterp<dirZ, InterpDir::PrimalToDual>();

            auto const& By_at_Ez = ct_state.getBy_at_Ez();
            auto const& Bz_at_Ey = ct_state.getBz_at_Ey();

            auto const EzBy
                = project_on_face(zEdgeToXFace, [&](auto const& i) { return Ez(i) * By_at_Ez(i); });
            auto const EyBz
                = project_on_face(yEdgeToXFace, [&](auto const& i) { return Ey(i) * Bz_at_Ey(i); });

            F_Etot += EyBz - EzBy;
        }
        else if constexpr (direction == Direction::Y)
        {
            constexpr auto zEdgeToYFace
                = GridLayout::implT::template directionalInterp<dirX, InterpDir::PrimalToDual>();
            constexpr auto xEdgeToYFace
                = GridLayout::implT::template directionalInterp<dirZ, InterpDir::PrimalToDual>();

            auto const& Bx_at_Ez = ct_state.getBx_at_Ez();
            auto const& Bz_at_Ex = ct_state.getBz_at_Ex();

            auto const EzBx
                = project_on_face(zEdgeToYFace, [&](auto const& i) { return Ez(i) * Bx_at_Ez(i); });
            auto const ExBz
                = project_on_face(xEdgeToYFace, [&](auto const& i) { return Ex(i) * Bz_at_Ex(i); });

            F_Etot += EzBx - ExBz;
        }
        else if constexpr (direction == Direction::Z)
        {
            constexpr auto xEdgeToZFace
                = GridLayout::implT::template directionalInterp<dirY, InterpDir::PrimalToDual>();
            constexpr auto yEdgeToZFace
                = GridLayout::implT::template directionalInterp<dirX, InterpDir::PrimalToDual>();

            auto const& By_at_Ex = ct_state.getBy_at_Ex();
            auto const& Bx_at_Ey = ct_state.getBx_at_Ey();

            auto const ExBy
                = project_on_face(xEdgeToZFace, [&](auto const& i) { return Ex(i) * By_at_Ex(i); });
            auto const EyBx
                = project_on_face(yEdgeToZFace, [&](auto const& i) { return Ey(i) * Bx_at_Ey(i); });

            F_Etot += ExBy - EyBx;
        }
    }


    GridLayout layout_;
    Equations equations_;
    RiemannSolver_t riemann_;
};

} // namespace PHARE::core

#endif
