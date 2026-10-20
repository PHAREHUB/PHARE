#ifndef PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP
#define PHARE_CORE_NUMERICS_GODUNOV_FLUXES_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/reconstructions/reconstructor.hpp"

#include "initializer/data_provider.hpp"

#include <array>
#include <limits>
#include <vector>
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

            riemann_fluxes_<direction>(ct_state, state, fluxes);

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
    // Face states come from a cell pass: each cell reconstructs the values at both of its faces
    // from one stencil (Reconstruction::cell_faces), and the face pass reads its left cell's
    // upper value and its right cell's lower value. Fields that are not cell-centred (transverse
    // B, J) are projected to cell centres once beforehand.
    template<auto direction>
    void riemann_fluxes_(auto& ct_state, auto const& state, auto& fluxes)
    {
        using implT                = GridLayout::implT;
        constexpr auto d           = static_cast<std::size_t>(direction);
        constexpr int ng           = Reconstruction_t::nghosts;
        constexpr std::size_t nCt  = Hall ? 5 : 2; // projected to cell centres: Bt0 Bt1 (Jx Jy Jz)
        constexpr std::size_t nRec = 5 + nCt;      // rho Vx Vy Vz P, then the projected ones

        constexpr auto transverse = []() {
            if constexpr (direction == Direction::X)
                return std::array{Direction::Y, Direction::Z};
            else if constexpr (direction == Direction::Y)
                return std::array{Direction::X, Direction::Z};
            else
                return std::array{Direction::X, Direction::Y};
        }();

        auto const& faces = fluxes.template expose_centering<direction>();
        std::array<int, dimension> lo, hi;
        for (std::size_t i = 0; i < dimension; ++i)
        {
            auto const [start, end] = layout_.physicalStartToEnd(faces, static_cast<Direction>(i));
            lo[i]                   = static_cast<int>(start) - (i == d ? 0 : ng);
            hi[i]                   = static_cast<int>(end) + (i == d ? 0 : ng);
        }
        auto cellLo = lo;
        cellLo[d] -= 1; // cells on both sides of the faces
        auto projLo = cellLo, projHi = hi;
        projLo[d] -= ng - 1; // cell_faces stencil reach
        projHi[d] += ng - 1;

        auto const shape    = state.rho.shape();
        std::size_t const N = product(shape, std::size_t{1});
        static thread_local std::vector<double> buffer;
        if (buffer.size() < N * (nCt + 2 * nRec))
            buffer.resize(N * (nCt + 2 * nRec));
        auto view = [&](std::size_t slot) {
            return NdArrayView<dimension, double>{buffer.data() + slot * N, shape};
        };
        auto centred = [&](std::size_t v) { return view(v); };
        auto lower   = [&](std::size_t v) { return view(nCt + v); };
        auto upper   = [&](std::size_t v) { return view(nCt + nRec + v); };

        for_box_(projLo, projHi, [&](MeshIndex<dimension> const& c) {
            centred(0)(c)
                = GridLayout::template project<Reconstructor_t::template projection<transverse[0]>>(
                    state.B(static_cast<Component>(transverse[0])), c);
            centred(1)(c)
                = GridLayout::template project<Reconstructor_t::template projection<transverse[1]>>(
                    state.B(static_cast<Component>(transverse[1])), c);
            if constexpr (Hall)
            {
                centred(2)(c) = GridLayout::template project<implT::edgeXToCellCenter>(
                    state.J(Component::X), c);
                centred(3)(c) = GridLayout::template project<implT::edgeYToCellCenter>(
                    state.J(Component::Y), c);
                centred(4)(c) = GridLayout::template project<implT::edgeZToCellCenter>(
                    state.J(Component::Z), c);
            }
        });

        for_box_(cellLo, hi, [&](MeshIndex<dimension> const& c) {
            auto store = [&](std::size_t v, auto const& F) {
                auto const [l, u] = Reconstruction_t::template cell_faces<direction>(F, c);
                lower(v)(c)       = l;
                upper(v)(c)       = u;
            };
            store(0, state.rho);
            store(1, state.V(Component::X));
            store(2, state.V(Component::Y));
            store(3, state.V(Component::Z));
            store(4, state.P);
            for (std::size_t v = 0; v < nCt; ++v)
                store(5 + v, centred(v));
        });

        auto const& Bn = state.B(static_cast<Component>(direction));
        for_box_(lo, hi, [&](MeshIndex<dimension> const& f) {
            auto const cl = GridLayout::template previous<direction>(f);
            auto face     = [&](auto&& side, MeshIndex<dimension> const& c) {
                PerIndexVector<double> B;
                B(direction)     = Bn(f);
                B(transverse[0]) = side(5)(c);
                B(transverse[1]) = side(6)(c);
                return PerIndex{side(0)(c), {side(1)(c), side(2)(c), side(3)(c)}, B, side(4)(c)};
            };
            auto uL  = face(upper, cl);
            auto uR  = face(lower, f);
            auto&& u = std::forward_as_tuple(uL, uR);

            if constexpr (Hall)
            {
                PerIndexVector<double> jL{upper(7)(cl), upper(8)(cl), upper(9)(cl)};
                PerIndexVector<double> jR{lower(7)(f), lower(8)(f), lower(9)(f)};
                auto const& j = std::forward_as_tuple(jL, jR);

                auto const& [fL, fR] = for_N<2, for_N_R_mode::make_tuple>([&](auto i) {
                    return equations_.template compute<direction>(std::get<i>(u), std::get<i>(j));
                });

                fluxes.template get_dir<direction>(f)
                    = riemann_.template solve<direction>(uL, uR, fL, fR, jL, jR);

                ct_state.template save<direction>(riemann_.vt, riemann_.jt, riemann_.rhot,
                                                  riemann_.uct_coefs, f);
            }
            else
            {
                auto const& [fL, fR] = for_N<2, for_N_R_mode::make_tuple>(
                    [&](auto i) { return equations_.template compute<direction>(std::get<i>(u)); });

                fluxes.template get_dir<direction>(f)
                    = riemann_.template solve<direction>(uL, uR, fL, fR);

                ct_state.template save<direction>(riemann_.vt, riemann_.uct_coefs, f);
            }
        });
    }

    static void for_box_(auto const& lo, auto const& hi, auto&& fn)
    {
        for (int i = lo[0]; i <= hi[0]; ++i)
            if constexpr (dimension == 1)
                fn(make_index(i));
            else
                for (int j = lo[1]; j <= hi[1]; ++j)
                    if constexpr (dimension == 2)
                        fn(make_index(i, j));
                    else
                        for (int k = lo[2]; k <= hi[2]; ++k)
                            fn(make_index(i, j, k));
    }

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
