#ifndef CORE_NUMERICS_RECONSTRUCTION_WENOZ_HPP
#define CORE_NUMERICS_RECONSTRUCTION_WENOZ_HPP

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include "core/numerics/reconstructions/stencil.hpp"
#include <utility>

namespace PHARE::core
{
template<typename GridLayout, typename SlopeLimiter = void>
class WENOZReconstruction
{
public:
    static constexpr auto nghosts = 3;

    using GridLayout_t = GridLayout;

    template<auto direction, typename Field>
    static auto reconstruct(Field const& F, MeshIndex<Field::dimension> index)
    {
        auto u_3
            = F(GridLayout::template previous<direction>(GridLayout::template previous<direction>(
                GridLayout::template previous<direction>(index))));
        auto u_2 = F(GridLayout::template previous<direction>(
            GridLayout::template previous<direction>(index)));
        auto u_1 = F(GridLayout::template previous<direction>(index));
        auto u   = F(index);
        auto u1  = F(GridLayout::template next<direction>(index));
        auto u2
            = F(GridLayout::template next<direction>(GridLayout::template next<direction>(index)));

        return std::make_pair(recons_wenoz_L_(u_3, u_2, u_1, u, u1),
                              recons_wenoz_R_(u_2, u_1, u, u1, u2));
    }

    template<auto direction, auto src, auto dst, typename Field>
    static auto center_reconstruct(Field const& U, MeshIndex<Field::dimension> index)
    {
        auto u_3 = GridLayout::template project<src, dst>(
            U, GridLayout::template previous<direction>(GridLayout::template previous<direction>(
                   GridLayout::template previous<direction>(index))));
        auto u_2 = GridLayout::template project<src, dst>(
            U, GridLayout::template previous<direction>(
                   GridLayout::template previous<direction>(index)));
        auto u_1 = GridLayout::template project<src, dst>(
            U, GridLayout::template previous<direction>(index));
        auto u  = GridLayout::template project<src, dst>(U, index);
        auto u1 = GridLayout::template project<src, dst>(
            U, GridLayout::template next<direction>(index));
        auto u2 = GridLayout::template project<src, dst>(
            U, GridLayout::template next<direction>(GridLayout::template next<direction>(index)));

        return std::make_pair(recons_wenoz_L_(u_3, u_2, u_1, u, u1),
                              recons_wenoz_R_(u_2, u_1, u, u1, u2));
    }

    // values at the lower and upper faces of cell `index`, from the one stencil centred on it
    template<auto direction, typename Field>
    static auto cell_faces(Field const& F, MeshIndex<Field::dimension> index)
    {
        return faces_of(centred_stencil<GridLayout, nghosts - 1, direction>(F, index));
    }

    // lower and upper face values of the centre cell of a 2 * (nghosts - 1) + 1 stencil
    static auto faces_of(auto const& u)
    {
        // both faces share the scaled weights x_k, only the linear weights d_k differ; one
        // divide serves both normalisations: 1/sR = sL / (sR * sL)
        auto const [x0, x1, x2] = wenoz_scaled_(u[0], u[1], u[2], u[3], u[4]);

        auto const sR = dR0 * x0 + dR1 * x1 + dR2 * x2;
        auto const sL = dL0 * x0 + dL1 * x1 + dL2 * x2;
        auto const nR = recons_wenoz_R_from_(u[0], u[1], u[2], u[3], u[4],
                                             std::make_tuple(dR0 * x0, dR1 * x1, dR2 * x2));
        auto const nL = recons_wenoz_L_from_(u[0], u[1], u[2], u[3], u[4],
                                             std::make_tuple(dL0 * x0, dL1 * x1, dL2 * x2));

        auto const inv = 1. / (sR * sL);
        return std::make_pair(nR * (sL * inv), nL * (sR * inv));
    }

private:
    static constexpr auto dL0 = 1. / 10.;
    static constexpr auto dL1 = 3. / 5.;
    static constexpr auto dL2 = 3. / 10.;
    static constexpr auto dR0 = 3. / 10.;
    static constexpr auto dR1 = 3. / 5.;
    static constexpr auto dR2 = 1. / 10.;

    static auto recons_wenoz_L_(auto const ull, auto const ul, auto const u, auto const ur,
                                auto const urr)
    {
        auto const [x0, x1, x2] = wenoz_scaled_(ull, ul, u, ur, urr);
        return recons_wenoz_L_from_(ull, ul, u, ur, urr, normalized_(x0, x1, x2, dL0, dL1, dL2));
    }

    static auto recons_wenoz_R_(auto const ull, auto const ul, auto const u, auto const ur,
                                auto const urr)
    {
        auto const [x0, x1, x2] = wenoz_scaled_(ull, ul, u, ur, urr);
        return recons_wenoz_R_from_(ull, ul, u, ur, urr, normalized_(x0, x1, x2, dR0, dR1, dR2));
    }

    static auto recons_wenoz_L_from_(auto const ull, auto const ul, auto const u, auto const ur,
                                     auto const urr, auto const& w)
    {
        auto const [wL0, wL1, wL2] = w;

        return wL0 * ((1. / 3.) * ull - (7. / 6.) * ul + (11. / 6.) * u)
               + wL1 * (-(1. / 6.) * ul + (5. / 6.) * u + (1. / 3.) * ur)
               + wL2 * ((1. / 3.) * u + (5. / 6.) * ur - (1. / 6.) * urr);
    }

    static auto recons_wenoz_R_from_(auto const ull, auto const ul, auto const u, auto const ur,
                                     auto const urr, auto const& w)
    {
        auto const [wR0, wR1, wR2] = w;

        return wR0 * ((1. / 3.) * u + (5. / 6.) * ul - (1. / 6.) * ull)
               + wR1 * (-(1. / 6.) * ur + (5. / 6.) * u + (1. / 3.) * ul)
               + wR2 * ((1. / 3.) * urr - (7. / 6.) * ur + (11. / 6.) * u);
    }

    // x_k = (1 + tau5 / (beta_k + eps)) * prod_j (beta_j + eps): the stencil-dependent part of
    // the WENO-Z weights, independent of which face is reconstructed. The common factor
    // prod_j (beta_j + eps) cancels in the normalisation and removes the divides. It stays
    // representable since beta_j + eps >= 1e-40 and beta_j ~ |du|^2 (|du| << 1e25).
    static auto wenoz_scaled_(auto const ull, auto const ul, auto const u, auto const ur,
                              auto const urr)
    {
        static constexpr auto eps = 1.e-40;

        auto const beta0 = (13. / 12.) * (ull - 2. * ul + u) * (ull - 2. * ul + u)
                           + (1. / 4.) * (ull - 4. * ul + 3. * u) * (ull - 4. * ul + 3. * u);
        auto const beta1 = (13. / 12.) * (ul - 2. * u + ur) * (ul - 2. * u + ur)
                           + (1. / 4.) * (ul - ur) * (ul - ur);
        auto const beta2 = (13. / 12.) * (u - 2. * ur + urr) * (u - 2. * ur + urr)
                           + (1. / 4.) * (3. * u - 4. * ur + urr) * (3. * u - 4. * ur + urr);

        auto const tau5 = std::abs(beta0 - beta2);

        auto const D0 = beta0 + eps;
        auto const D1 = beta1 + eps;
        auto const D2 = beta2 + eps;

        return std::make_tuple((D0 + tau5) * (D1 * D2), (D1 + tau5) * (D0 * D2),
                               (D2 + tau5) * (D0 * D1));
    }

    static auto normalized_(auto const r0, auto const r1, auto const r2, auto const d0,
                            auto const d1, auto const d2)
    {
        auto const alpha0 = d0 * r0;
        auto const alpha1 = d1 * r1;
        auto const alpha2 = d2 * r2;

        auto const inv_sum = 1. / (alpha0 + alpha1 + alpha2);

        return std::make_tuple(alpha0 * inv_sum, alpha1 * inv_sum, alpha2 * inv_sum);
    }
};

} // namespace PHARE::core

#endif
