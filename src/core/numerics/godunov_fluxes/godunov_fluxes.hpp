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
#include "core/numerics/reconstructions/stencil.hpp"

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
        // cells per line along the last dimension; also the variable stride in a line buffer
        int const n = hi[dimension - 1] - cellLo[dimension - 1] + 1;
        // projected fields over the whole patch, then two line buffers of nRec lower and nRec
        // upper face values, so face states stay in cache between the cell and face passes
        static thread_local std::vector<double> buffer;
        if (buffer.size() < N * nCt + 4 * nRec * n)
            buffer.resize(N * nCt + 4 * nRec * n);
        double* const scratch = buffer.data(); // one thread_local lookup, not one per access
        std::array<double*, 2> lines{scratch + N * nCt, scratch + N * nCt + 2 * nRec * n};
        auto view = [&](std::size_t slot) {
            return NdArrayView<dimension, double>{scratch + slot * N, shape};
        };
        auto centred = [&](std::size_t v) { return view(v); }; // then nRec lower, nRec upper

        for_each_index<dimension>(projLo, projHi, [&](MeshIndex<dimension> const& c) {
            centred(0)(c)
                = GridLayout::template project<Reconstructor_t::template projection<transverse[0]>>(
                    state.B[static_cast<Component>(transverse[0])], c);
            centred(1)(c)
                = GridLayout::template project<Reconstructor_t::template projection<transverse[1]>>(
                    state.B[static_cast<Component>(transverse[1])], c);
            if constexpr (Hall)
            {
                centred(2)(c) = GridLayout::template project<implT::edgeXToCellCenter>(
                    state.J[Component::X], c);
                centred(3)(c) = GridLayout::template project<implT::edgeYToCellCenter>(
                    state.J[Component::Y], c);
                centred(4)(c) = GridLayout::template project<implT::edgeZToCellCenter>(
                    state.J[Component::Z], c);
            }
        });

        // cell pass of one contiguous line (last dimension) starting at cell c, into a line buffer
        std::array<double const*, nRec> src{state.rho.data(), state.V[Component::X].data(),
                                            state.V[Component::Y].data(),
                                            state.V[Component::Z].data(), state.P.data()};
        for (std::size_t v = 0; v < nCt; ++v)
            src[5 + v] = scratch + v * N;
        std::ptrdiff_t stride = 1; // between neighbours along direction, all arrays share shape
        for (std::size_t i = d + 1; i < dimension; ++i)
            stride *= shape[i];
        auto cell_line = [&](MeshIndex<dimension> const& c, double* line) {
            auto const o = NdArrayViewer<dimension>::idx(shape, c);
            for (std::size_t v = 0; v < nRec; ++v)
                cell_faces_line<Reconstruction_t>(src[v] + o, stride, n, line + v * n,
                                                  line + (nRec + v) * n);
        };

        // face pass of one contiguous line (last dimension) starting at face f0, from the upper
        // values of its left cells (up) and the lower values of its right cells (lw). Every array
        // is read and written through a pointer to the line start, so consecutive faces are
        // consecutive elements whatever each array's shape; the Riemann solver is copied so the
        // per-face state it keeps in members (vt, jt, rhot, uct_coefs) lives in registers.
        auto const& Bn = state.B[static_cast<Component>(direction)];
        auto const ct  = [&]() {
            if constexpr (direction == Direction::X)
                return std::forward_as_tuple(ct_state.vt_x, ct_state.aL_x, ct_state.aR_x,
                                              ct_state.dL_x, ct_state.dR_x);
            else if constexpr (direction == Direction::Y)
                return std::forward_as_tuple(ct_state.vt_y, ct_state.aL_y, ct_state.aR_y,
                                              ct_state.dL_y, ct_state.dR_y);
            else
                return std::forward_as_tuple(ct_state.vt_z, ct_state.aL_z, ct_state.aR_z,
                                              ct_state.dL_z, ct_state.dR_z);
        }();
        auto faceHi           = hi;
        faceHi[dimension - 1] = lo[dimension - 1];
        int const nf          = hi[dimension - 1] - lo[dimension - 1] + 1;

        auto face_line = [&](MeshIndex<dimension> const& f0, double const* up, double const* lw) {
            double const* bn = &Bn(f0);

            auto F0      = fluxes.template get_dir<direction>(f0);
            double* fRho = &F0.rho;
            double* fVx  = &F0.V.x;
            double* fVy  = &F0.V.y;
            double* fVz  = &F0.V.z;
            double* fE   = &F0.P;

            auto& [vT, aL, aR, dL, dR] = ct;
            double* vtx                = &vT[Component::X](f0);
            double* vty                = &vT[Component::Y](f0);
            double* vtz                = &vT[Component::Z](f0);
            double* aLp                = &aL(f0);
            double* aRp                = &aR(f0);
            double* dLp                = &dL(f0);
            double* dRp                = &dR(f0);

            auto riemann = riemann_;

            auto state_at = [&](double const* p, int k) PHARE_LAMBDA_ALWAYS_INLINE {
                auto const B = [&]() -> PerIndexVector<double> {
                    if constexpr (direction == Direction::X)
                        return {bn[k], p[5 * n + k], p[6 * n + k]};
                    else if constexpr (direction == Direction::Y)
                        return {p[5 * n + k], bn[k], p[6 * n + k]};
                    else
                        return {p[5 * n + k], p[6 * n + k], bn[k]};
                }();
                return PerIndex{p[k], {p[n + k], p[2 * n + k], p[3 * n + k]}, B, p[4 * n + k]};
            };
            auto store = [&](int k, auto const& F) PHARE_LAMBDA_ALWAYS_INLINE {
                fRho[k] = F.rho;
                fVx[k]  = F.V.x;
                fVy[k]  = F.V.y;
                fVz[k]  = F.V.z;
                fE[k]   = F.P;
                vtx[k]  = riemann.vt.x;
                vty[k]  = riemann.vt.y;
                vtz[k]  = riemann.vt.z;
                aLp[k]  = riemann.uct_coefs[0];
                aRp[k]  = riemann.uct_coefs[1];
                dLp[k]  = riemann.uct_coefs[2];
                dRp[k]  = riemann.uct_coefs[3];
            };

            if constexpr (Hall)
            {
                auto& jT     = ct_state.template getJt<direction>();
                double* jtx  = &jT[Component::X](f0);
                double* jty  = &jT[Component::Y](f0);
                double* jtz  = &jT[Component::Z](f0);
                double* rhot = &ct_state.template getRhot<direction>()(f0);

                PHARE_VECTORIZE_IVDEP
                for (int k = 0; k < nf; ++k)
                {
                    auto uL = state_at(up, k);
                    auto uR = state_at(lw, k);
                    PerIndexVector<double> jL{up[7 * n + k], up[8 * n + k], up[9 * n + k]};
                    PerIndexVector<double> jR{lw[7 * n + k], lw[8 * n + k], lw[9 * n + k]};

                    auto const fL = equations_.template compute<direction>(uL, jL);
                    auto const fR = equations_.template compute<direction>(uR, jR);

                    store(k, riemann.template solve<direction>(uL, uR, fL, fR, jL, jR));
                    jtx[k]  = riemann.jt.x;
                    jty[k]  = riemann.jt.y;
                    jtz[k]  = riemann.jt.z;
                    rhot[k] = riemann.rhot;
                }
            }
            else
            {
                PHARE_VECTORIZE_IVDEP
                for (int k = 0; k < nf; ++k)
                {
                    auto uL = state_at(up, k);
                    auto uR = state_at(lw, k);

                    auto const fL = equations_.template compute<direction>(uL);
                    auto const fR = equations_.template compute<direction>(uR);

                    store(k, riemann.template solve<direction>(uL, uR, fL, fR));
                }
            }
        };

        if constexpr (d == dimension - 1)
        {
            // faces and their two cells share a line: face k sits between cells k and k + 1
            for_each_index<dimension>(lo, faceHi, [&](MeshIndex<dimension> const& f0) {
                cell_line(GridLayout::template previous<direction>(f0), lines[0]);
                face_line(f0, lines[0] + nRec * n, lines[0] + 1);
            });
        }
        else
        {
            // the lines along direction are walked in order, so a face line reads the upper
            // values of the previous cell line and the lower values of the current one
            auto outerHi = faceHi;
            outerHi[d]   = lo[d];
            for_each_index<dimension>(lo, outerHi, [&](MeshIndex<dimension> const& f) {
                auto c = GridLayout::template previous<direction>(f);
                cell_line(c, lines[0]);
                for (int i = lo[d]; i <= hi[d]; ++i)
                {
                    c[d] = i;
                    std::swap(lines[0], lines[1]);
                    cell_line(c, lines[0]);
                    face_line(c, lines[1] + nRec * n, lines[0]);
                }
            });
        }
    }

    template<auto direction>
    auto transverse_on_face_(auto const& E, auto const& B, MeshIndex<dimension> idx) const
    {
        using implT        = GridLayout::implT;
        auto constexpr nan = std::numeric_limits<double>::quiet_NaN();

        auto const& Ex = E[Component::X];
        auto const& Ey = E[Component::Y];
        auto const& Ez = E[Component::Z];
        auto const& Bx = B[Component::X];
        auto const& By = B[Component::Y];
        auto const& Bz = B[Component::Z];

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
