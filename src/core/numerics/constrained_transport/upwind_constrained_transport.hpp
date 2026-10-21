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

#include <array>
#include <cmath>
#include <vector>

namespace PHARE::core
{
using UpwindConstrainedTransportInfo = OhmInfo;

template<typename GridLayout, template<typename> typename Reconstruction, bool Hall>
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
        // scratch holds one component's slots; slot numbers as in ExEq_ EyEq_ EzEq_
        constexpr std::size_t perDir = Hall ? 4 : 2; // slots per reconstruction direction
        constexpr std::size_t nSlots = dimension == 1 ? 1 : 2 * perDir;
        auto const N
            = std::max({product(Ex.shape(), std::size_t{1}), product(Ey.shape(), std::size_t{1}),
                        product(Ez.shape(), std::size_t{1})});
        static thread_local std::vector<double> buffer;
        if (buffer.size() < 2 * nSlots * N)
            buffer.resize(2 * nSlots * N);

        double* const scratch = buffer.data(); // one thread_local lookup, not one per access

        auto slots = [&](auto const& Ec) {
            return [&](std::size_t slot, MeshIndex<dimension> const& idx)
                       PHARE_LAMBDA_ALWAYS_INLINE {
                           auto const s = 2 * slot * N;
                           return std::make_pair(
                               NdArrayView<dimension, double>{scratch + s, Ec.shape()}(idx),
                               NdArrayView<dimension, double>{scratch + s + N, Ec.shape()}(idx));
                       };
        };
        // slots sharing a direction go through one cell loop, so their independent stencils
        // overlap instead of each running as its own latency chain
        auto const transverse = [&]<auto D, auto Dt>(auto const& Ec, std::size_t first,
                                                     auto const& Vt, auto const& Bc) {
            if constexpr (Hall)
                face_states_<D>(Ec, scratch, N, first, Vt, Bc,
                                ct_state.template getJt<Dt>()[static_cast<Component>(D)],
                                ct_state.template getRhot<Dt>());
            else
                face_states_<D>(Ec, scratch, N, first, Vt, Bc);
        };

        if constexpr (dimension == 2)
            face_states_<Direction::Y>(Ex, scratch, N, 0, B[Component::Z]);
        else if constexpr (dimension == 3)
        {
            transverse.template operator()<Direction::Y, Direction::Z>(
                Ex, 0, ct_state.vt_z[Component::Y], B[Component::Z]);
            transverse.template operator()<Direction::Z, Direction::Y>(
                Ex, perDir, ct_state.vt_y[Component::Z], B[Component::Y]);
        }
        if constexpr (dimension > 1)
        {
            auto const exFaces = slots(Ex);
            layout_.evalOnBox(Ex,
                              [&](auto&... args) { ExEq_(ct_state, Ex, B, {args...}, exFaces); });
        }

        if constexpr (dimension <= 2)
            face_states_<Direction::X>(Ey, scratch, N, 0, B[Component::Z]);
        else
        {
            transverse.template operator()<Direction::X, Direction::Z>(
                Ey, 0, ct_state.vt_z[Component::X], B[Component::Z]);
            transverse.template operator()<Direction::Z, Direction::X>(
                Ey, perDir, ct_state.vt_x[Component::Z], B[Component::X]);
        }
        auto const eyFaces = slots(Ey);
        layout_.evalOnBox(Ey, [&](auto&... args) { EyEq_(ct_state, Ey, B, {args...}, eyFaces); });

        if constexpr (dimension == 1)
            face_states_<Direction::X>(Ez, scratch, N, 0, B[Component::Y]);
        else
        {
            transverse.template operator()<Direction::Y, Direction::X>(
                Ez, 0, ct_state.vt_x[Component::Y], B[Component::X]);
            transverse.template operator()<Direction::X, Direction::Y>(
                Ez, perDir, ct_state.vt_y[Component::X], B[Component::Y]);
        }
        auto const ezFaces = slots(Ez);
        layout_.evalOnBox(Ez, [&](auto&... args) { EzEq_(ct_state, Ez, B, {args...}, ezFaces); });

        if (is_dissipative_)
        {
            auto const& Ediss = dissipative_electric_state.E();
            for (auto c : {Component::X, Component::Y, Component::Z})
                layout_.evalOnBox(E(c), [&](auto&... args) { E(c)(args...) += Ediss(c)(args...); });
        }
    }

private:
    // For every index of E's physical box, the reconstructed values of each F on the edge from
    // the cell below it (slot's first half, L) and above it (second half, R) along D. Each cell's
    // stencil is evaluated once and gives both of its faces, instead of once per adjacent edge.
    template<auto D, typename... Fields>
    void face_states_(auto const& E, double* scratch, std::size_t N, std::size_t first,
                      Fields const&... F) const
    {
        std::array<int, dimension> lo, hi;
        for (std::size_t i = 0; i < dimension; ++i)
        {
            auto const [start, end] = layout_.physicalStartToEnd(E, static_cast<Direction>(i));
            lo[i]                   = static_cast<int>(start);
            hi[i]                   = static_cast<int>(end);
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
                        scratch + 2 * slot * N + N + atC, scratch + 2 * slot * N + atAbove);
                    ++slot;
                }(),
                ...);
        });
    }

    void ExEq_(auto& ct_state, auto& Ex, auto const& B, MeshIndex<dimension> idx,
               auto const& faces) const
    {
        if constexpr (dimension == 2)
        {
            auto [BzL, BzR] = faces(0, idx);

            auto FL = BzL * ct_state.vt_y[Component::Y](idx)
                      - B[Component::Y](idx) * ct_state.vt_y[Component::Z](idx);
            auto FR = BzR * ct_state.vt_y[Component::Y](idx)
                      - B[Component::Y](idx) * ct_state.vt_y[Component::Z](idx);

            Ex(idx) = -(ct_state.aL_y(idx) * FL + ct_state.aR_y(idx) * FR)
                      + (ct_state.dR_y(idx) * BzR - ct_state.dL_y(idx) * BzL);

            if constexpr (Hall)
            {
                auto invRho  = 1.0 / ct_state.template getRhot<Direction::Y>()(idx);
                auto JxB_x_L = ct_state.template getJt<Direction::Y>()[Component::Y](idx) * BzL
                               - ct_state.template getJt<Direction::Y>()[Component::Z](idx)
                                     * B(Component::Y)(idx);
                auto JxB_x_R = ct_state.template getJt<Direction::Y>()[Component::Y](idx) * BzR
                               - ct_state.template getJt<Direction::Y>()[Component::Z](idx)
                                     * B(Component::Y)(idx);

                auto HallL = -JxB_x_L * invRho;
                auto HallR = -JxB_x_R * invRho;

                auto F_Bz_y = ct_state.aL_y(idx) * HallL + ct_state.aR_y(idx) * HallR;

                Ex(idx) += -F_Bz_y;
            }
        }
        else if constexpr (dimension == 3)
        {
            auto aS = 0.5
                      * (ct_state.aL_y(idx)
                         + ct_state.aL_y(layout_.template previous<Direction::Z>(idx)));
            auto aN = 0.5
                      * (ct_state.aR_y(idx)
                         + ct_state.aR_y(layout_.template previous<Direction::Z>(idx)));
            auto aB = 0.5
                      * (ct_state.aL_z(idx)
                         + ct_state.aL_z(layout_.template previous<Direction::Y>(idx)));
            auto aT = 0.5
                      * (ct_state.aR_z(idx)
                         + ct_state.aR_z(layout_.template previous<Direction::Y>(idx)));
            auto dS = 0.5
                      * (ct_state.dL_y(idx)
                         + ct_state.dL_y(layout_.template previous<Direction::Z>(idx)));
            auto dN = 0.5
                      * (ct_state.dR_y(idx)
                         + ct_state.dR_y(layout_.template previous<Direction::Z>(idx)));
            auto dB = 0.5
                      * (ct_state.dL_z(idx)
                         + ct_state.dL_z(layout_.template previous<Direction::Y>(idx)));
            auto dT = 0.5
                      * (ct_state.dR_z(idx)
                         + ct_state.dR_z(layout_.template previous<Direction::Y>(idx)));

            constexpr std::size_t z0 = Hall ? 4 : 2; // first slot reconstructed along Z
            auto [vyS, vyN]          = faces(0, idx);
            auto [vzB, vzT]          = faces(z0, idx);
            auto [BzS, BzN]          = faces(1, idx);
            auto [ByB, ByT]          = faces(z0 + 1, idx);

            Ex(idx) = (aB * vzB * ByB + aT * vzT * ByT) - (aS * vyS * BzS + aN * vyN * BzN)
                      - (dT * ByT - dB * ByB) + (dN * BzN - dS * BzS);

            if constexpr (Hall)
            {
                auto [jyS, jyN]   = faces(2, idx);
                auto [jzB, jzT]   = faces(z0 + 2, idx);
                auto [rhoS, rhoN] = faces(3, idx);
                auto [rhoB, rhoT] = faces(z0 + 3, idx);

                Ex(idx) += -(aB * jzB * ByB / rhoB + aT * jzT * ByT / rhoT)
                           + (aS * jyS * BzS / rhoS + aN * jyN * BzN / rhoN);
            }
        }
    }

    void EyEq_(auto& ct_state, auto& Ey, auto const& B, MeshIndex<dimension> idx,
               auto const& faces) const
    {
        if constexpr (dimension <= 2)
        {
            auto [BzL, BzR] = faces(0, idx);

            auto FL = BzL * ct_state.vt_x[Component::X](idx)
                      - B[Component::X](idx) * ct_state.vt_x[Component::Z](idx);
            auto FR = BzR * ct_state.vt_x[Component::X](idx)
                      - B[Component::X](idx) * ct_state.vt_x[Component::Z](idx);

            Ey(idx) = (ct_state.aL_x(idx) * FL + ct_state.aR_x(idx) * FR)
                      - (ct_state.dR_x(idx) * BzR - ct_state.dL_x(idx) * BzL);

            if constexpr (Hall)
            {
                auto invRho  = 1.0 / ct_state.template getRhot<Direction::X>()(idx);
                auto JxB_y_L = ct_state.template getJt<Direction::X>()[Component::Z](idx)
                                   * B(Component::X)(idx)
                               - ct_state.template getJt<Direction::X>()[Component::X](idx) * BzL;
                auto JxB_y_R = ct_state.template getJt<Direction::X>()[Component::Z](idx)
                                   * B(Component::X)(idx)
                               - ct_state.template getJt<Direction::X>()[Component::X](idx) * BzR;

                auto HallL = JxB_y_L * invRho;
                auto HallR = JxB_y_R * invRho;

                auto F_Bz_x = ct_state.aL_x(idx) * HallL + ct_state.aR_x(idx) * HallR;

                Ey(idx) += F_Bz_x;
            }
        }
        else if constexpr (dimension == 3)
        {
            auto aW = 0.5
                      * (ct_state.aL_x(idx)
                         + ct_state.aL_x(layout_.template previous<Direction::Z>(idx)));
            auto aE = 0.5
                      * (ct_state.aR_x(idx)
                         + ct_state.aR_x(layout_.template previous<Direction::Z>(idx)));
            auto aB = 0.5
                      * (ct_state.aL_z(idx)
                         + ct_state.aL_z(layout_.template previous<Direction::X>(idx)));
            auto aT = 0.5
                      * (ct_state.aR_z(idx)
                         + ct_state.aR_z(layout_.template previous<Direction::X>(idx)));
            auto dW = 0.5
                      * (ct_state.dL_x(idx)
                         + ct_state.dL_x(layout_.template previous<Direction::Z>(idx)));
            auto dE = 0.5
                      * (ct_state.dR_x(idx)
                         + ct_state.dR_x(layout_.template previous<Direction::Z>(idx)));
            auto dB = 0.5
                      * (ct_state.dL_z(idx)
                         + ct_state.dL_z(layout_.template previous<Direction::X>(idx)));
            auto dT = 0.5
                      * (ct_state.dR_z(idx)
                         + ct_state.dR_z(layout_.template previous<Direction::X>(idx)));

            constexpr std::size_t z0 = Hall ? 4 : 2; // first slot reconstructed along Z
            auto [vxW, vxE]          = faces(0, idx);
            auto [vzB, vzT]          = faces(z0, idx);
            auto [BzW, BzE]          = faces(1, idx);
            auto [BxB, BxT]          = faces(z0 + 1, idx);

            Ey(idx) = (aW * vxW * BzW + aE * vxE * BzE) - (aB * vzB * BxB + aT * vzT * BxT)
                      - (dE * BzE - dW * BzW) + (dT * BxT - dB * BxB);

            if constexpr (Hall)
            {
                auto [jxW, jxE]   = faces(2, idx);
                auto [jzB, jzT]   = faces(z0 + 2, idx);
                auto [rhoW, rhoE] = faces(3, idx);
                auto [rhoB, rhoT] = faces(z0 + 3, idx);
                Ey(idx) += -(aW * jxW * BzW / rhoW + aE * jxE * BzE / rhoE)
                           + (aB * jzB * BxB / rhoB + aT * jzT * BxT / rhoT);
            }
        }
    }

    void EzEq_(auto& ct_state, auto& Ez, auto const& B, MeshIndex<dimension> idx,
               auto const& faces) const
    {
        if constexpr (dimension == 1)
        {
            auto [ByL, ByR] = faces(0, idx);

            auto FL = ByL * ct_state.vt_x[Component::X](idx)
                      - B[Component::X](idx) * ct_state.vt_x[Component::Y](idx);
            auto FR = ByR * ct_state.vt_x[Component::X](idx)
                      - B[Component::X](idx) * ct_state.vt_x[Component::Y](idx);

            Ez(idx) = -(ct_state.aL_x(idx) * FL + ct_state.aR_x(idx) * FR)
                      + (ct_state.dR_x(idx) * ByR - ct_state.dL_x(idx) * ByL);

            if constexpr (Hall)
            {
                auto invRho  = 1.0 / ct_state.template getRhot<Direction::X>()(idx);
                auto JxB_z_L = ct_state.template getJt<Direction::X>()[Component::X](idx) * ByL
                               - ct_state.template getJt<Direction::X>()[Component::Y](idx)
                                     * B(Component::X)(idx);
                auto JxB_z_R = ct_state.template getJt<Direction::X>()[Component::X](idx) * ByR
                               - ct_state.template getJt<Direction::X>()[Component::Y](idx)
                                     * B(Component::X)(idx);

                auto HallL = -JxB_z_L * invRho;
                auto HallR = -JxB_z_R * invRho;

                auto F_By_x = ct_state.aL_x(idx) * HallL + ct_state.aR_x(idx) * HallR;

                Ez(idx) += -F_By_x;
            }
        }
        else if constexpr (dimension >= 2)
        {
            auto aW = 0.5
                      * (ct_state.aL_x(idx)
                         + ct_state.aL_x(layout_.template previous<Direction::Y>(idx)));
            auto aE = 0.5
                      * (ct_state.aR_x(idx)
                         + ct_state.aR_x(layout_.template previous<Direction::Y>(idx)));
            auto aS = 0.5
                      * (ct_state.aL_y(idx)
                         + ct_state.aL_y(layout_.template previous<Direction::X>(idx)));
            auto aN = 0.5
                      * (ct_state.aR_y(idx)
                         + ct_state.aR_y(layout_.template previous<Direction::X>(idx)));
            auto dW = 0.5
                      * (ct_state.dL_x(idx)
                         + ct_state.dL_x(layout_.template previous<Direction::Y>(idx)));
            auto dE = 0.5
                      * (ct_state.dR_x(idx)
                         + ct_state.dR_x(layout_.template previous<Direction::Y>(idx)));
            auto dS = 0.5
                      * (ct_state.dL_y(idx)
                         + ct_state.dL_y(layout_.template previous<Direction::X>(idx)));
            auto dN = 0.5
                      * (ct_state.dR_y(idx)
                         + ct_state.dR_y(layout_.template previous<Direction::X>(idx)));

            constexpr std::size_t x0 = Hall ? 4 : 2; // first slot reconstructed along X
            auto [vyS, vyN]          = faces(0, idx);
            auto [vxW, vxE]          = faces(x0, idx);
            auto [BxS, BxN]          = faces(1, idx);
            auto [ByW, ByE]          = faces(x0 + 1, idx);

            Ez(idx) = -(aW * vxW * ByW + aE * vxE * ByE) + (aS * vyS * BxS + aN * vyN * BxN)
                      + (dE * ByE - dW * ByW) - (dN * BxN - dS * BxS);

            if constexpr (Hall)
            {
                auto [jyS, jyN]   = faces(2, idx);
                auto [jxW, jxE]   = faces(x0 + 2, idx);
                auto [rhoS, rhoN] = faces(3, idx);
                auto [rhoW, rhoE] = faces(x0 + 3, idx);

                Ez(idx) += (aW * jxW * ByW / rhoW + aE * jxE * ByE / rhoE)
                           - (aS * jyS * BxS / rhoS + aN * jyN * BxN / rhoN);
            }
        }
    }

    GridLayout layout_;
    bool const is_dissipative_;
};
} // namespace PHARE::core

#endif
