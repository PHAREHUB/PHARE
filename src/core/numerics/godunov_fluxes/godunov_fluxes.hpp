#ifndef PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP
#define PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/reconstructions/reconstructor.hpp"

#include "initializer/data_provider.hpp"

#include <limits>
#include <utility>
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

template<auto direction, size_t dim>
auto getGrow(int const nghosts)
{
    Point<std::uint32_t, dim> p{};

    auto dir = static_cast<size_t>(direction);

    for (size_t i = 0; i < dim; ++i)
        if (i != dir)
            p[i] = nghosts;

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

    constexpr static auto Hall = Equations::hall;

    static_assert(GridLayout::options.field_ghost_width >= Reconstruction_t::nghosts + 1,
                  "MHD ghost width too small for the reconstruction stencil plus ampere's layer");

    explicit Godunov(GodunovInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
        , is_resistive_{info.isResistive()}
        , is_hyper_resistive_{info.isHyperResistive()}
        , equations_{gamma}
        , riemann_{gamma}
    {
    }

    template<typename State, typename Fluxes>
    void operator()(auto& ct_state, auto const& dissipative_electric_state, State& state,
                    Fluxes& fluxes)
    {
        constexpr auto directions = getDirections<dimension>();

        constexpr auto num_directions = std::tuple_size_v<std::decay_t<decltype(directions)>>;

        for_N<num_directions>([&](auto i) {
            constexpr Direction direction = std::get<i>(directions);

            layout_.evalOnBiggerBox(
                fluxes.template expose_centering<direction>(),
                getGrow<direction, dimension>(Reconstruction_t::nghosts), [&](auto&... indices) {
                    auto&& [uL, uR]
                        = Reconstructor_t::template reconstruct<direction>(state, {indices...});
                    auto&& u = std::forward_as_tuple(uL, uR);

                    if constexpr (Hall)
                    {
                        auto const& [jL, jR] = Reconstructor_t::template center_reconstruct<
                            direction, GridLayout::implT::edgeXToCellCenter,
                            GridLayout::implT::edgeYToCellCenter,
                            GridLayout::implT::edgeZToCellCenter>(state.J, {indices...});

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
                    else
                    {
                        auto const& [fL, fR] = for_N<2, for_N_R_mode::make_tuple>([&](auto i) {
                            return equations_.template compute<direction>(std::get<i>(u));
                        });

                        fluxes.template get_dir<direction>({indices...})
                            = riemann_.template solve<direction>(uL, uR, fL, fR);

                        ct_state.template save<direction>(riemann_.vt, riemann_.uct_coefs,
                                                          {indices...});
                    }
                });

            if (is_resistive_ || is_hyper_resistive_)
                layout_.evalOnBox(fluxes.template expose_centering<direction>(),
                                  [&](auto&... indices) {
                                      auto F = fluxes.template get_dir<direction>({indices...});
                                      auto const [Et, Bt] = transverse_on_face_<direction>(
                                          dissipative_electric_state.E(), state.B, {indices...});
                                      equations_.template dissipative_contributions<direction>(
                                          Et, Bt, F.B, F.Etot());
                                  });
        });
    }

    bool isResistive() const { return is_resistive_; }
    bool isHyperResistive() const { return is_hyper_resistive_; }

private:
    template<auto direction>
    auto transverse_on_face_(auto const& E, auto const& B, MeshIndex<dimension> idx) const
    {
        using implT        = GridLayout::implT;
        auto constexpr nan = std::numeric_limits<double>::quiet_NaN();

        auto const& Ex = E(Component::X);
        auto const& Ey = E(Component::Y);
        auto const& Ez = E(Component::Z);
        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        if constexpr (direction == Direction::X)
            return std::make_pair(
                PerIndexVector<double>{nan,
                                       GridLayout::template project<implT::edgeYToFaceX>(Ey, idx),
                                       GridLayout::template project<implT::edgeZToFaceX>(Ez, idx)},
                PerIndexVector<double>{nan, GridLayout::template project<implT::ByToFaceX>(By, idx),
                                       GridLayout::template project<implT::BzToFaceX>(Bz, idx)});
        else if constexpr (direction == Direction::Y)
            return std::make_pair(
                PerIndexVector<double>{GridLayout::template project<implT::edgeXToFaceY>(Ex, idx),
                                       nan,
                                       GridLayout::template project<implT::edgeZToFaceY>(Ez, idx)},
                PerIndexVector<double>{GridLayout::template project<implT::BxToFaceY>(Bx, idx), nan,
                                       GridLayout::template project<implT::BzToFaceY>(Bz, idx)});
        else
            return std::make_pair(
                PerIndexVector<double>{GridLayout::template project<implT::edgeXToFaceZ>(Ex, idx),
                                       GridLayout::template project<implT::edgeYToFaceZ>(Ey, idx),
                                       nan},
                PerIndexVector<double>{GridLayout::template project<implT::BxToFaceZ>(Bx, idx),
                                       GridLayout::template project<implT::ByToFaceZ>(By, idx),
                                       nan});
    }

    GridLayout layout_;
    bool const is_resistive_;
    bool const is_hyper_resistive_;
    Equations equations_;
    RiemannSolver_t riemann_;
};

} // namespace PHARE::core

#endif
