#ifndef PHARE_UPWIND_CONSTRAINED_TRANSPORT_HPP
#define PHARE_UPWIND_CONSTRAINED_TRANSPORT_HPP


#include "core/def.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/index/index.hpp"

#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/utilities/types.hpp"
#include "core/numerics/reconstructions/stencil.hpp"
#include "core/numerics/riemann_solvers/uct_coefficients.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace PHARE::core
{
using UpwindConstrainedTransportInfo = OhmInfo;

template<typename GridLayout, template<typename> typename Reconstruction, bool Hall,
         UCTCoefficients uct = UCTCoefficients{}>
class UpwindConstrainedTransport : UpwindConstrainedTransportInfo
{
    using Super                     = UpwindConstrainedTransportInfo;
    using Reconstruction_t          = Reconstruction<GridLayout>;
    constexpr static auto dimension = GridLayout::dimension;
    using Super::hyper_mode;

public:
    using Info_t = Super;

    UpwindConstrainedTransport(UpwindConstrainedTransportInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
        , is_dissipative_{info.isResistive() || info.isHyperResistive()}
    {
    }

    void operator()(auto& ct_state, auto const& dissipative_electric_state, auto& mhd_state) const
    {
        auto& E       = mhd_state.E;
        auto const& B = mhd_state.B;

        auto& Ex = E[Component::X];
        auto& Ey = E[Component::Y];
        auto& Ez = E[Component::Z];

        // edge states from a cell pass (see face_states_), one E component at a time so the
        // scratch holds one component's slots; slot numbers as in ExLine_ EyLine_ EzLine_
        constexpr std::size_t perDir = Hall ? 4 : 2; // slots per reconstruction direction
        constexpr std::size_t nSlots = dimension == 1 ? 1 : 2 * perDir;
        static thread_local std::vector<double> buffer;

        // Each component goes through the scratch a slab of x-planes at a time: on big patches
        // the whole scratch outgrows the cache and would make a round trip to memory between
        // the cell pass writing it and the edge pass reading it back. A slab of nx edge planes
        // keeps nx + 2 planes of every slot (its X stencils reach one plane either side).
        constexpr std::size_t slabBytes = std::size_t{4} << 20;
        auto const forEachSlab          = [&](auto const& Ec, auto&& states, auto&& edges) {
            auto const shape  = Ec.shape();
            std::size_t plane = 1;
            for (std::size_t i = 1; i < dimension; ++i)
                plane *= shape[i];
            auto const [xStart, xEnd] = layout_.physicalStartToEnd(Ec, Direction::X);
            int const start = static_cast<int>(xStart), end = static_cast<int>(xEnd);
            int nx = end - start + 1;
            if constexpr (dimension > 1)
                nx = std::clamp(static_cast<int>(slabBytes / (2 * nSlots * plane * sizeof(double))),
                                         1, nx);
            std::size_t const Ns = static_cast<std::size_t>(nx + 2) * plane; // per half slot
            if (buffer.size() < 2 * nSlots * Ns)
                buffer.resize(2 * nSlots * Ns);
            double* const scratch = buffer.data(); // one thread_local lookup, not one per access

            for (int x0 = start; x0 <= end; x0 += nx)
            {
                Slab const slab{x0, std::min(x0 + nx - 1, end),
                                static_cast<std::size_t>(x0 - 1) * plane, Ns, scratch};
                states(slab);
                edges(slab);
            }
        };

        // the line of Ec's edges from `first` in the slab, and its first edge in slot 0
        auto combine = [&](auto const& Ec, Slab const& slab, auto&& line) {
            std::array<int, dimension> lo, hi;
            for (std::size_t i = 0; i < dimension; ++i)
            {
                auto const [start, end] = layout_.physicalStartToEnd(Ec, static_cast<Direction>(i));
                lo[i]                   = static_cast<int>(start);
                hi[i]                   = static_cast<int>(end);
            }
            auto const n = static_cast<std::uint32_t>(hi[dimension - 1] - lo[dimension - 1] + 1);
            if constexpr (dimension > 1)
            {
                lo[0] = slab.x0;
                hi[0] = slab.x1;
            }
            hi[dimension - 1] = lo[dimension - 1];
            for_each_index<dimension>(lo, hi, [&](MeshIndex<dimension> const& first) {
                line(first, n, slab.at(NdArrayViewer<dimension>::idx(Ec.shape(), first)), slab.N);
            });
        };

        // slots sharing a direction go through one cell loop, so their independent stencils
        // overlap instead of each running as its own latency chain
        auto const transverse
            = [&]<auto D, auto Dt>(auto const& Ec, Slab const& slab, std::size_t first,
                                   auto const& Vt, auto const& Bc) {
                  if constexpr (Hall)
                      face_states_<D>(Ec, slab, first, Vt, Bc,
                                      ct_state.template getJt<Dt>()[static_cast<Component>(D)],
                                      ct_state.template getRhot<Dt>());
                  else
                      face_states_<D>(Ec, slab, first, Vt, Bc);
              };

        if constexpr (dimension > 1)
            forEachSlab(
                Ex,
                [&](Slab const& slab) {
                    if constexpr (dimension == 2)
                        face_states_<Direction::Y>(Ex, slab, 0, B[Component::Z]);
                    else
                    {
                        transverse.template operator()<Direction::Y, Direction::Z>(
                            Ex, slab, 0, ct_state.vt_z[Component::Y], B[Component::Z]);
                        transverse.template operator()<Direction::Z, Direction::Y>(
                            Ex, slab, perDir, ct_state.vt_y[Component::Z], B[Component::Y]);
                    }
                },
                [&](Slab const& slab) {
                    combine(Ex, slab, [&](auto const& first, auto n, double const* f, auto N) {
                        ExLine_(ct_state, Ex, B, first, n, f, N);
                    });
                });

        forEachSlab(
            Ey,
            [&](Slab const& slab) {
                if constexpr (dimension <= 2)
                    face_states_<Direction::X>(Ey, slab, 0, B[Component::Z]);
                else
                {
                    transverse.template operator()<Direction::X, Direction::Z>(
                        Ey, slab, 0, ct_state.vt_z[Component::X], B[Component::Z]);
                    transverse.template operator()<Direction::Z, Direction::X>(
                        Ey, slab, perDir, ct_state.vt_x[Component::Z], B[Component::X]);
                }
            },
            [&](Slab const& slab) {
                combine(Ey, slab, [&](auto const& first, auto n, double const* f, auto N) {
                    EyLine_(ct_state, Ey, B, first, n, f, N);
                });
            });

        forEachSlab(
            Ez,
            [&](Slab const& slab) {
                if constexpr (dimension == 1)
                    face_states_<Direction::X>(Ez, slab, 0, B[Component::Y]);
                else
                {
                    transverse.template operator()<Direction::Y, Direction::X>(
                        Ez, slab, 0, ct_state.vt_x[Component::Y], B[Component::X]);
                    transverse.template operator()<Direction::X, Direction::Y>(
                        Ez, slab, perDir, ct_state.vt_y[Component::X], B[Component::Y]);
                }
            },
            [&](Slab const& slab) {
                combine(Ez, slab, [&](auto const& first, auto n, double const* f, auto N) {
                    EzLine_(ct_state, Ez, B, first, n, f, N);
                });
            });

        if (is_dissipative_)
        {
            auto const& Ediss = dissipative_electric_state.E();
            for_N<3>([&](auto i) {
                constexpr auto c = static_cast<Component>(i());
                auto& Ec         = E(c);
                auto const& Edc  = Ediss(c);
                layout_.evalOnLines(Ec,
                                    [&](MeshIndex<dimension> const& first, std::uint32_t const n) {
                                        double* ec        = &Ec(first);
                                        double const* edc = &Edc(first);
                                        PHARE_VECTORIZE_IVDEP
                                        for (std::uint32_t k = 0; k < n; ++k)
                                            ec[k] += edc[k];
                                    });
            });
        }
    }

private:
    // edge planes x0..x1 of one E component, and its scratch: a half slot holds N values from
    // the plane before x0 on, so that the edge at flat index i of E is at(i)
    struct Slab
    {
        int x0, x1;
        std::size_t base, N;
        double* scratch;

        double* at(std::size_t i) const { return scratch + (i - base); }
    };

    // For every index of the slab's part of E's physical box, the reconstructed values of each F
    // on the edge from the cell below it (slot's first half, L) and above it (second half, R)
    // along D. Each cell's stencil is evaluated once and gives both of its faces, instead of once
    // per adjacent edge. Kept out of line: inlined, GCC vectorises it with masked loops, which
    // run slower.
    template<auto D, typename... Fields>
    [[gnu::noinline]] void face_states_(auto const& E, Slab const& slab, std::size_t first,
                                        Fields const&... F) const
    {
        std::array<int, dimension> lo, hi;
        for (std::size_t i = 0; i < dimension; ++i)
        {
            auto const [start, end] = layout_.physicalStartToEnd(E, static_cast<Direction>(i));
            lo[i]                   = static_cast<int>(start);
            hi[i]                   = static_cast<int>(end);
        }
        if constexpr (dimension > 1)
        {
            lo[0] = slab.x0;
            hi[0] = slab.x1;
        }
        lo[static_cast<std::size_t>(D)] -= 1;

        // one contiguous line (last dimension) at a time per field
        auto const shape      = E.shape();
        auto lineHi           = hi;
        lineHi[dimension - 1] = lo[dimension - 1];
        int const n           = hi[dimension - 1] - lo[dimension - 1] + 1;
        for_each_index<dimension>(lo, lineHi, [&](MeshIndex<dimension> const& c) {
            auto const atC = NdArrayViewer<dimension>::idx(shape, c);
            auto const atAbove
                = NdArrayViewer<dimension>::idx(shape, GridLayout::template next<D>(c));
            std::size_t slot = first;
            (
                [&] {
                    auto const fShape     = F.shape();
                    std::ptrdiff_t stride = 1;
                    for (std::size_t i = static_cast<std::size_t>(D) + 1; i < dimension; ++i)
                        stride *= fShape[i];
                    cell_faces_line<Reconstruction_t>(
                        F.data() + NdArrayViewer<dimension>::idx(fShape, c), stride, n,
                        slab.at(atC) + 2 * slot * slab.N + slab.N,
                        slab.at(atAbove) + 2 * slot * slab.N);
                    ++slot;
                }(),
                ...);
        });
    }

    // ExLine_ EyLine_ EzLine_: the E component on the line of n edges from `first`,
    // f pointing at the line's first edge in slot 0 of the scratch (slot j: states
    // from the cell below the edge at f + 2 j N, from the cell above N further).
    // Every field is read through a pointer to its own line, so that the loop
    // vectorises.
    static double const* at_(auto const& F, MeshIndex<dimension> const& index) { return &F(index); }

    // 1/2 at every edge: the advective coefficients of a solver that does not store them
    struct Half
    {
        PHARE_ALWAYS_INLINE double operator[](std::uint32_t) const { return 0.5; }
    };

    // the line of an advective coefficient (aL, aR) and of the right diffusive one (dR), as the
    // Riemann solver stores them (see UCTCoefficients)
    static auto adv_(auto const& a, MeshIndex<dimension> const& index)
    {
        if constexpr (uct.halfAdvective)
            return Half{};
        else
            return at_(a, index);
    }
    static double const* diffR_(auto const& dL, auto const& dR, MeshIndex<dimension> const& index)
    {
        if constexpr (uct.equalDiffusive)
            return at_(dL, index);
        else
            return at_(dR, index);
    }

    void ExLine_(auto& ct_state, auto& Ex, auto const& B, MeshIndex<dimension> const& first,
                 std::uint32_t const n, double const* f, std::size_t const N) const
    {
        auto const L = [&](std::size_t slot) { return f + 2 * slot * N; };
        auto const R = [&](std::size_t slot) { return f + 2 * slot * N + N; };
        double* ex   = &Ex(first);

        if constexpr (dimension == 2)
        {
            auto const *BzL = L(0), *BzR = R(0);
            auto const* vy   = at_(ct_state.vt_y[Component::Y], first);
            auto const* vz   = at_(ct_state.vt_y[Component::Z], first);
            auto const* By   = at_(B[Component::Y], first);
            auto const aL    = adv_(ct_state.aL_y, first);
            auto const aR    = adv_(ct_state.aR_y, first);
            auto const* dL   = at_(ct_state.dL_y, first);
            auto const* dR   = diffR_(ct_state.dL_y, ct_state.dR_y, first);
            double const *jy = nullptr, *jz = nullptr, *rho = nullptr; // Hall only
            if constexpr (Hall)
            {
                auto const& jt = ct_state.template getJt<Direction::Y>();
                jy             = at_(jt[Component::Y], first);
                jz             = at_(jt[Component::Z], first);
                rho            = at_(ct_state.template getRhot<Direction::Y>(), first);
            }

            PHARE_VECTORIZE_IVDEP
            for (std::uint32_t k = 0; k < n; ++k)
            {
                auto FL = BzL[k] * vy[k] - By[k] * vz[k];
                auto FR = BzR[k] * vy[k] - By[k] * vz[k];

                double e = -(aL[k] * FL + aR[k] * FR) + (dR[k] * BzR[k] - dL[k] * BzL[k]);

                if constexpr (Hall)
                {
                    auto invRho  = 1.0 / rho[k];
                    auto JxB_x_L = jy[k] * BzL[k] - jz[k] * By[k];
                    auto JxB_x_R = jy[k] * BzR[k] - jz[k] * By[k];

                    auto HallL = -JxB_x_L * invRho;
                    auto HallR = -JxB_x_R * invRho;

                    auto F_Bz_y = aL[k] * HallL + aR[k] * HallR;

                    e += -F_Bz_y;
                }
                ex[k] = e;
            }
        }
        else if constexpr (dimension == 3)
        {
            auto const pY = layout_.template previous<Direction::Y>(first);
            auto const pZ = layout_.template previous<Direction::Z>(first);

            auto const aLy = adv_(ct_state.aL_y, first), aLy_ = adv_(ct_state.aL_y, pZ);
            auto const aRy = adv_(ct_state.aR_y, first), aRy_ = adv_(ct_state.aR_y, pZ);
            auto const aLz = adv_(ct_state.aL_z, first), aLz_ = adv_(ct_state.aL_z, pY);
            auto const aRz = adv_(ct_state.aR_z, first), aRz_ = adv_(ct_state.aR_z, pY);
            auto const *dLy = at_(ct_state.dL_y, first), *dLy_ = at_(ct_state.dL_y, pZ);
            auto const *dRy  = diffR_(ct_state.dL_y, ct_state.dR_y, first),
                       *dRy_ = diffR_(ct_state.dL_y, ct_state.dR_y, pZ);
            auto const *dLz = at_(ct_state.dL_z, first), *dLz_ = at_(ct_state.dL_z, pY);
            auto const *dRz  = diffR_(ct_state.dL_z, ct_state.dR_z, first),
                       *dRz_ = diffR_(ct_state.dL_z, ct_state.dR_z, pY);

            constexpr std::size_t z0 = Hall ? 4 : 2; // first slot reconstructed along Z

            PHARE_VECTORIZE_IVDEP
            for (std::uint32_t k = 0; k < n; ++k)
            {
                auto aS = 0.5 * (aLy[k] + aLy_[k]);
                auto aN = 0.5 * (aRy[k] + aRy_[k]);
                auto aB = 0.5 * (aLz[k] + aLz_[k]);
                auto aT = 0.5 * (aRz[k] + aRz_[k]);
                auto dS = 0.5 * (dLy[k] + dLy_[k]);
                auto dN = 0.5 * (dRy[k] + dRy_[k]);
                auto dB = 0.5 * (dLz[k] + dLz_[k]);
                auto dT = 0.5 * (dRz[k] + dRz_[k]);

                auto vyS = L(0)[k], vyN = R(0)[k];
                auto vzB = L(z0)[k], vzT = R(z0)[k];
                auto BzS = L(1)[k], BzN = R(1)[k];
                auto ByB = L(z0 + 1)[k], ByT = R(z0 + 1)[k];

                double e = (aB * vzB * ByB + aT * vzT * ByT) - (aS * vyS * BzS + aN * vyN * BzN)
                           - (dT * ByT - dB * ByB) + (dN * BzN - dS * BzS);

                if constexpr (Hall)
                {
                    auto jyS = L(2)[k], jyN = R(2)[k];
                    auto jzB = L(z0 + 2)[k], jzT = R(z0 + 2)[k];
                    auto rhoS = L(3)[k], rhoN = R(3)[k];
                    auto rhoB = L(z0 + 3)[k], rhoT = R(z0 + 3)[k];

                    e += -(aB * jzB * ByB / rhoB + aT * jzT * ByT / rhoT)
                         + (aS * jyS * BzS / rhoS + aN * jyN * BzN / rhoN);
                }
                ex[k] = e;
            }
        }
    }

    void EyLine_(auto& ct_state, auto& Ey, auto const& B, MeshIndex<dimension> const& first,
                 std::uint32_t const n, double const* f, std::size_t const N) const
    {
        auto const L = [&](std::size_t slot) { return f + 2 * slot * N; };
        auto const R = [&](std::size_t slot) { return f + 2 * slot * N + N; };
        double* ey   = &Ey(first);

        if constexpr (dimension <= 2)
        {
            auto const *BzL = L(0), *BzR = R(0);
            auto const* vx   = at_(ct_state.vt_x[Component::X], first);
            auto const* vz   = at_(ct_state.vt_x[Component::Z], first);
            auto const* Bx   = at_(B[Component::X], first);
            auto const aL    = adv_(ct_state.aL_x, first);
            auto const aR    = adv_(ct_state.aR_x, first);
            auto const* dL   = at_(ct_state.dL_x, first);
            auto const* dR   = diffR_(ct_state.dL_x, ct_state.dR_x, first);
            double const *jx = nullptr, *jz = nullptr, *rho = nullptr; // Hall only
            if constexpr (Hall)
            {
                auto const& jt = ct_state.template getJt<Direction::X>();
                jx             = at_(jt[Component::X], first);
                jz             = at_(jt[Component::Z], first);
                rho            = at_(ct_state.template getRhot<Direction::X>(), first);
            }

            PHARE_VECTORIZE_IVDEP
            for (std::uint32_t k = 0; k < n; ++k)
            {
                auto FL = BzL[k] * vx[k] - Bx[k] * vz[k];
                auto FR = BzR[k] * vx[k] - Bx[k] * vz[k];

                double e = (aL[k] * FL + aR[k] * FR) - (dR[k] * BzR[k] - dL[k] * BzL[k]);

                if constexpr (Hall)
                {
                    auto invRho  = 1.0 / rho[k];
                    auto JxB_y_L = jz[k] * Bx[k] - jx[k] * BzL[k];
                    auto JxB_y_R = jz[k] * Bx[k] - jx[k] * BzR[k];

                    auto HallL = JxB_y_L * invRho;
                    auto HallR = JxB_y_R * invRho;

                    auto F_Bz_x = aL[k] * HallL + aR[k] * HallR;

                    e += F_Bz_x;
                }
                ey[k] = e;
            }
        }
        else if constexpr (dimension == 3)
        {
            auto const pX = layout_.template previous<Direction::X>(first);
            auto const pZ = layout_.template previous<Direction::Z>(first);

            auto const aLx = adv_(ct_state.aL_x, first), aLx_ = adv_(ct_state.aL_x, pZ);
            auto const aRx = adv_(ct_state.aR_x, first), aRx_ = adv_(ct_state.aR_x, pZ);
            auto const aLz = adv_(ct_state.aL_z, first), aLz_ = adv_(ct_state.aL_z, pX);
            auto const aRz = adv_(ct_state.aR_z, first), aRz_ = adv_(ct_state.aR_z, pX);
            auto const *dLx = at_(ct_state.dL_x, first), *dLx_ = at_(ct_state.dL_x, pZ);
            auto const *dRx  = diffR_(ct_state.dL_x, ct_state.dR_x, first),
                       *dRx_ = diffR_(ct_state.dL_x, ct_state.dR_x, pZ);
            auto const *dLz = at_(ct_state.dL_z, first), *dLz_ = at_(ct_state.dL_z, pX);
            auto const *dRz  = diffR_(ct_state.dL_z, ct_state.dR_z, first),
                       *dRz_ = diffR_(ct_state.dL_z, ct_state.dR_z, pX);

            constexpr std::size_t z0 = Hall ? 4 : 2; // first slot reconstructed along Z

            PHARE_VECTORIZE_IVDEP
            for (std::uint32_t k = 0; k < n; ++k)
            {
                auto aW = 0.5 * (aLx[k] + aLx_[k]);
                auto aE = 0.5 * (aRx[k] + aRx_[k]);
                auto aB = 0.5 * (aLz[k] + aLz_[k]);
                auto aT = 0.5 * (aRz[k] + aRz_[k]);
                auto dW = 0.5 * (dLx[k] + dLx_[k]);
                auto dE = 0.5 * (dRx[k] + dRx_[k]);
                auto dB = 0.5 * (dLz[k] + dLz_[k]);
                auto dT = 0.5 * (dRz[k] + dRz_[k]);

                auto vxW = L(0)[k], vxE = R(0)[k];
                auto vzB = L(z0)[k], vzT = R(z0)[k];
                auto BzW = L(1)[k], BzE = R(1)[k];
                auto BxB = L(z0 + 1)[k], BxT = R(z0 + 1)[k];

                double e = (aW * vxW * BzW + aE * vxE * BzE) - (aB * vzB * BxB + aT * vzT * BxT)
                           - (dE * BzE - dW * BzW) + (dT * BxT - dB * BxB);

                if constexpr (Hall)
                {
                    auto jxW = L(2)[k], jxE = R(2)[k];
                    auto jzB = L(z0 + 2)[k], jzT = R(z0 + 2)[k];
                    auto rhoW = L(3)[k], rhoE = R(3)[k];
                    auto rhoB = L(z0 + 3)[k], rhoT = R(z0 + 3)[k];

                    e += -(aW * jxW * BzW / rhoW + aE * jxE * BzE / rhoE)
                         + (aB * jzB * BxB / rhoB + aT * jzT * BxT / rhoT);
                }
                ey[k] = e;
            }
        }
    }

    void EzLine_(auto& ct_state, auto& Ez, auto const& B, MeshIndex<dimension> const& first,
                 std::uint32_t const n, double const* f, std::size_t const N) const
    {
        auto const L = [&](std::size_t slot) { return f + 2 * slot * N; };
        auto const R = [&](std::size_t slot) { return f + 2 * slot * N + N; };
        double* ez   = &Ez(first);

        if constexpr (dimension == 1)
        {
            auto const *ByL = L(0), *ByR = R(0);
            auto const* vx   = at_(ct_state.vt_x[Component::X], first);
            auto const* vy   = at_(ct_state.vt_x[Component::Y], first);
            auto const* Bx   = at_(B[Component::X], first);
            auto const aL    = adv_(ct_state.aL_x, first);
            auto const aR    = adv_(ct_state.aR_x, first);
            auto const* dL   = at_(ct_state.dL_x, first);
            auto const* dR   = diffR_(ct_state.dL_x, ct_state.dR_x, first);
            double const *jx = nullptr, *jy = nullptr, *rho = nullptr; // Hall only
            if constexpr (Hall)
            {
                auto const& jt = ct_state.template getJt<Direction::X>();
                jx             = at_(jt[Component::X], first);
                jy             = at_(jt[Component::Y], first);
                rho            = at_(ct_state.template getRhot<Direction::X>(), first);
            }

            PHARE_VECTORIZE_IVDEP
            for (std::uint32_t k = 0; k < n; ++k)
            {
                auto FL = ByL[k] * vx[k] - Bx[k] * vy[k];
                auto FR = ByR[k] * vx[k] - Bx[k] * vy[k];

                double e = -(aL[k] * FL + aR[k] * FR) + (dR[k] * ByR[k] - dL[k] * ByL[k]);

                if constexpr (Hall)
                {
                    auto invRho  = 1.0 / rho[k];
                    auto JxB_z_L = jx[k] * ByL[k] - jy[k] * Bx[k];
                    auto JxB_z_R = jx[k] * ByR[k] - jy[k] * Bx[k];

                    auto HallL = -JxB_z_L * invRho;
                    auto HallR = -JxB_z_R * invRho;

                    auto F_By_x = aL[k] * HallL + aR[k] * HallR;

                    e += -F_By_x;
                }
                ez[k] = e;
            }
        }
        else if constexpr (dimension >= 2)
        {
            auto const pX = layout_.template previous<Direction::X>(first);
            auto const pY = layout_.template previous<Direction::Y>(first);

            auto const aLx = adv_(ct_state.aL_x, first), aLx_ = adv_(ct_state.aL_x, pY);
            auto const aRx = adv_(ct_state.aR_x, first), aRx_ = adv_(ct_state.aR_x, pY);
            auto const aLy = adv_(ct_state.aL_y, first), aLy_ = adv_(ct_state.aL_y, pX);
            auto const aRy = adv_(ct_state.aR_y, first), aRy_ = adv_(ct_state.aR_y, pX);
            auto const *dLx = at_(ct_state.dL_x, first), *dLx_ = at_(ct_state.dL_x, pY);
            auto const *dRx  = diffR_(ct_state.dL_x, ct_state.dR_x, first),
                       *dRx_ = diffR_(ct_state.dL_x, ct_state.dR_x, pY);
            auto const *dLy = at_(ct_state.dL_y, first), *dLy_ = at_(ct_state.dL_y, pX);
            auto const *dRy  = diffR_(ct_state.dL_y, ct_state.dR_y, first),
                       *dRy_ = diffR_(ct_state.dL_y, ct_state.dR_y, pX);

            constexpr std::size_t x0 = Hall ? 4 : 2; // first slot reconstructed along X

            PHARE_VECTORIZE_IVDEP
            for (std::uint32_t k = 0; k < n; ++k)
            {
                auto aW = 0.5 * (aLx[k] + aLx_[k]);
                auto aE = 0.5 * (aRx[k] + aRx_[k]);
                auto aS = 0.5 * (aLy[k] + aLy_[k]);
                auto aN = 0.5 * (aRy[k] + aRy_[k]);
                auto dW = 0.5 * (dLx[k] + dLx_[k]);
                auto dE = 0.5 * (dRx[k] + dRx_[k]);
                auto dS = 0.5 * (dLy[k] + dLy_[k]);
                auto dN = 0.5 * (dRy[k] + dRy_[k]);

                auto vyS = L(0)[k], vyN = R(0)[k];
                auto vxW = L(x0)[k], vxE = R(x0)[k];
                auto BxS = L(1)[k], BxN = R(1)[k];
                auto ByW = L(x0 + 1)[k], ByE = R(x0 + 1)[k];

                double e = -(aW * vxW * ByW + aE * vxE * ByE) + (aS * vyS * BxS + aN * vyN * BxN)
                           + (dE * ByE - dW * ByW) - (dN * BxN - dS * BxS);

                if constexpr (Hall)
                {
                    auto jyS = L(2)[k], jyN = R(2)[k];
                    auto jxW = L(x0 + 2)[k], jxE = R(x0 + 2)[k];
                    auto rhoS = L(3)[k], rhoN = R(3)[k];
                    auto rhoW = L(x0 + 3)[k], rhoE = R(x0 + 3)[k];

                    e += (aW * jxW * ByW / rhoW + aE * jxE * ByE / rhoE)
                         - (aS * jyS * BxS / rhoS + aN * jyN * BxN / rhoN);
                }
                ez[k] = e;
            }
        }
    }
    GridLayout layout_;
    bool const is_dissipative_;
};
} // namespace PHARE::core

#endif
