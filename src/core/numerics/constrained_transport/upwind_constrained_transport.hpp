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

        if constexpr (dimension == 2)
        {
            // edge states from a cell pass (see face_states_), slot numbers as in ExEq_ EyEq_ EzEq_
            constexpr std::size_t nSlots = Hall ? 10 : 6;
            auto const N = product(Ez.shape(), std::size_t{1}); // primal-primal: largest 2D field
            static thread_local std::vector<double> buffer;
            if (buffer.size() < 2 * nSlots * N)
                buffer.resize(2 * nSlots * N);

            double* const scratch = buffer.data(); // one thread_local lookup, not one per access

            auto slots = [&](auto const& E, std::size_t first) {
                return [&, first](std::size_t slot, MeshIndex<dimension> const& idx)
                           PHARE_LAMBDA_ALWAYS_INLINE {
                               auto const s = 2 * (first + slot) * N;
                               return std::make_pair(
                                   NdArrayView<dimension, double>{scratch + s, E.shape()}(idx),
                                   NdArrayView<dimension, double>{scratch + s + N, E.shape()}(idx));
                           };
            };
            // slots sharing a direction and an edge box go through one cell loop, so their
            // independent stencils overlap instead of each running as its own latency chain
            face_states_<Direction::Y>(Ex, scratch, N, 0, B[Component::Z]);
            face_states_<Direction::X>(Ey, scratch, N, 1, B[Component::Z]);
            if constexpr (Hall)
            {
                face_states_<Direction::Y>(Ez, scratch, N, 2, ct_state.vt_x[Component::Y],
                                           B[Component::X],
                                           ct_state.template getJt<Direction::X>()[Component::Y],
                                           ct_state.template getRhot<Direction::X>());
                face_states_<Direction::X>(Ez, scratch, N, 6, ct_state.vt_y[Component::X],
                                           B[Component::Y],
                                           ct_state.template getJt<Direction::Y>()[Component::X],
                                           ct_state.template getRhot<Direction::Y>());
            }
            else
            {
                face_states_<Direction::Y>(Ez, scratch, N, 2, ct_state.vt_x[Component::Y],
                                           B[Component::X]);
                face_states_<Direction::X>(Ez, scratch, N, 4, ct_state.vt_y[Component::X],
                                           B[Component::Y]);
            }

            auto const exFaces = slots(Ex, 0);
            auto const eyFaces = slots(Ey, 1);
            auto const ezFaces = slots(Ez, 2);
            layout_.evalOnBox(Ex,
                              [&](auto&... args) { ExEq_(ct_state, Ex, B, {args...}, exFaces); });
            layout_.evalOnBox(Ey,
                              [&](auto&... args) { EyEq_(ct_state, Ey, B, {args...}, eyFaces); });
            layout_.evalOnBox(Ez,
                              [&](auto&... args) { EzEq_(ct_state, Ez, B, {args...}, ezFaces); });
        }
        else
        {
            auto const none
                = [](std::size_t, MeshIndex<dimension> const&) { return std::make_pair(0., 0.); };
            layout_.evalOnBox(Ex, [&](auto&... args) { ExEq_(ct_state, Ex, B, {args...}, none); });
            layout_.evalOnBox(Ey, [&](auto&... args) { EyEq_(ct_state, Ey, B, {args...}, none); });
            layout_.evalOnBox(Ez, [&](auto&... args) { EzEq_(ct_state, Ez, B, {args...}, none); });
        }

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

            auto [vyS, vyN] = Reconstruction_t::template reconstruct<Direction::Y>(
                ct_state.vt_z[Component::Y], idx);
            auto [vzB, vzT] = Reconstruction_t::template reconstruct<Direction::Z>(
                ct_state.vt_y[Component::Z], idx);

            auto [BzS, BzN]
                = Reconstruction_t::template reconstruct<Direction::Y>(B[Component::Z], idx);
            auto [ByB, ByT]
                = Reconstruction_t::template reconstruct<Direction::Z>(B[Component::Y], idx);

            Ex(idx) = (aB * vzB * ByB + aT * vzT * ByT) - (aS * vyS * BzS + aN * vyN * BzN)
                      - (dT * ByT - dB * ByB) + (dN * BzN - dS * BzS);

            if constexpr (Hall)
            {
                auto [jyS, jyN] = Reconstruction_t::template reconstruct<Direction::Y>(
                    ct_state.template getJt<Direction::Z>()[Component::Y], idx);
                auto [jzB, jzT] = Reconstruction_t::template reconstruct<Direction::Z>(
                    ct_state.template getJt<Direction::Y>()[Component::Z], idx);

                auto [rhoS, rhoN] = Reconstruction_t::template reconstruct<Direction::Y>(
                    ct_state.template getRhot<Direction::Z>(), idx);
                auto [rhoB, rhoT] = Reconstruction_t::template reconstruct<Direction::Z>(
                    ct_state.template getRhot<Direction::Y>(), idx);

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
            auto [BzL, BzR] = [&]() {
                if constexpr (dimension == 2)
                    return faces(0, idx);
                else
                    return Reconstruction_t::template reconstruct<Direction::X>(B[Component::Z],
                                                                                idx);
            }();

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

            auto [vxW, vxE] = Reconstruction_t::template reconstruct<Direction::X>(
                ct_state.vt_z[Component::X], idx);
            auto [vzB, vzT] = Reconstruction_t::template reconstruct<Direction::Z>(
                ct_state.vt_x[Component::Z], idx);
            auto [BzW, BzE]
                = Reconstruction_t::template reconstruct<Direction::X>(B[Component::Z], idx);
            auto [BxB, BxT]
                = Reconstruction_t::template reconstruct<Direction::Z>(B[Component::X], idx);

            Ey(idx) = (aW * vxW * BzW + aE * vxE * BzE) - (aB * vzB * BxB + aT * vzT * BxT)
                      - (dE * BzE - dW * BzW) + (dT * BxT - dB * BxB);

            if constexpr (Hall)
            {
                auto [jxW, jxE] = Reconstruction_t::template reconstruct<Direction::X>(
                    ct_state.template getJt<Direction::Z>()[Component::X], idx);
                auto [jzB, jzT] = Reconstruction_t::template reconstruct<Direction::Z>(
                    ct_state.template getJt<Direction::X>()[Component::Z], idx);
                auto [rhoW, rhoE] = Reconstruction_t::template reconstruct<Direction::X>(
                    ct_state.template getRhot<Direction::Z>(), idx);
                auto [rhoB, rhoT] = Reconstruction_t::template reconstruct<Direction::Z>(
                    ct_state.template getRhot<Direction::X>(), idx);
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
            auto [ByL, ByR]
                = Reconstruction_t::template reconstruct<Direction::X>(B[Component::Y], idx);

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

            auto const rec
                = [&]<auto D>(std::size_t slot, auto const& F) PHARE_LAMBDA_ALWAYS_INLINE {
                      if constexpr (dimension == 2)
                          return faces(slot, idx);
                      else
                          return Reconstruction_t::template reconstruct<D>(F, idx);
                  };
            constexpr std::size_t x0 = Hall ? 4 : 2; // first Ez slot reconstructed along X
            auto [vyS, vyN] = rec.template operator()<Direction::Y>(0, ct_state.vt_x[Component::Y]);
            auto [vxW, vxE]
                = rec.template operator()<Direction::X>(x0, ct_state.vt_y[Component::X]);

            auto [BxS, BxN] = rec.template operator()<Direction::Y>(1, B[Component::X]);
            auto [ByW, ByE] = rec.template operator()<Direction::X>(x0 + 1, B[Component::Y]);

            Ez(idx) = -(aW * vxW * ByW + aE * vxE * ByE) + (aS * vyS * BxS + aN * vyN * BxN)
                      + (dE * ByE - dW * ByW) - (dN * BxN - dS * BxS);

            if constexpr (Hall)
            {
                auto [jyS, jyN] = rec.template operator()<Direction::Y>(
                    2, ct_state.template getJt<Direction::X>()[Component::Y]);
                auto [jxW, jxE] = rec.template operator()<Direction::X>(
                    x0 + 2, ct_state.template getJt<Direction::Y>()[Component::X]);

                auto [rhoS, rhoN] = rec.template operator()<Direction::Y>(
                    3, ct_state.template getRhot<Direction::X>());
                auto [rhoW, rhoE] = rec.template operator()<Direction::X>(
                    x0 + 3, ct_state.template getRhot<Direction::Y>());

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
