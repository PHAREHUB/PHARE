#ifndef CORE_NUMERICS_RECONSTRUCTION_WENOZ_HPP
#define CORE_NUMERICS_RECONSTRUCTION_WENOZ_HPP

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
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

    template<auto direction, typename Field>
    static auto center_reconstruct(Field const& U, MeshIndex<Field::dimension> index,
                                   auto projection)
    {
        auto u_3 = GridLayout::project(
            U,
            GridLayout::template previous<direction>(GridLayout::template previous<direction>(
                GridLayout::template previous<direction>(index))),
            projection);
        auto u_2 = GridLayout::project(U,
                                       GridLayout::template previous<direction>(
                                           GridLayout::template previous<direction>(index)),
                                       projection);
        auto u_1
            = GridLayout::project(U, GridLayout::template previous<direction>(index), projection);
        auto u  = GridLayout::project(U, index, projection);
        auto u1 = GridLayout::project(U, GridLayout::template next<direction>(index), projection);
        auto u2 = GridLayout::project(
            U, GridLayout::template next<direction>(GridLayout::template next<direction>(index)),
            projection);

        return std::make_pair(recons_wenoz_L_(u_3, u_2, u_1, u, u1),
                              recons_wenoz_R_(u_2, u_1, u, u1, u2));
    }

private:
    static auto recons_wenoz_L_(auto const ull, auto const ul, auto const u, auto const ur,
                                auto const urr)
    {
        static constexpr auto dL0 = 1. / 10.;
        static constexpr auto dL1 = 3. / 5.;
        static constexpr auto dL2 = 3. / 10.;

        auto const [wL0, wL1, wL2] = compute_wenoz_weights(ull, ul, u, ur, urr, dL0, dL1, dL2);

        return wL0 * ((1. / 3.) * ull - (7. / 6.) * ul + (11. / 6.) * u)
               + wL1 * (-(1. / 6.) * ul + (5. / 6.) * u + (1. / 3.) * ur)
               + wL2 * ((1. / 3.) * u + (5. / 6.) * ur - (1. / 6.) * urr);
    }

    static auto recons_wenoz_R_(auto const ull, auto const ul, auto const u, auto const ur,
                                auto const urr)
    {
        static constexpr auto dR0 = 3. / 10.;
        static constexpr auto dR1 = 3. / 5.;
        static constexpr auto dR2 = 1. / 10.;

        auto const [wR0, wR1, wR2] = compute_wenoz_weights(ull, ul, u, ur, urr, dR0, dR1, dR2);

        return wR0 * ((1. / 3.) * u + (5. / 6.) * ul - (1. / 6.) * ull)
               + wR1 * (-(1. / 6.) * ur + (5. / 6.) * u + (1. / 3.) * ul)
               + wR2 * ((1. / 3.) * urr - (7. / 6.) * ur + (11. / 6.) * u);
    }

    static auto compute_wenoz_weights(auto const ull, auto const ul, auto const u, auto const ur,
                                      auto const urr, auto const d0, auto const d1, auto const d2)
    {
        static constexpr auto eps = 1.e-40;

        auto const beta0 = (13. / 12.) * (ull - 2. * ul + u) * (ull - 2. * ul + u)
                           + (1. / 4.) * (ull - 4. * ul + 3. * u) * (ull - 4. * ul + 3. * u);
        auto const beta1 = (13. / 12.) * (ul - 2. * u + ur) * (ul - 2. * u + ur)
                           + (1. / 4.) * (ul - ur) * (ul - ur);
        auto const beta2 = (13. / 12.) * (u - 2. * ur + urr) * (u - 2. * ur + urr)
                           + (1. / 4.) * (3. * u - 4. * ur + urr) * (3. * u - 4. * ur + urr);

        auto const tau5 = std::abs(beta0 - beta2);

        auto const alpha0 = d0 * (1. + tau5 / (beta0 + eps));
        auto const alpha1 = d1 * (1. + tau5 / (beta1 + eps));
        auto const alpha2 = d2 * (1. + tau5 / (beta2 + eps));

        auto const sum_alpha = alpha0 + alpha1 + alpha2;

        return std::make_tuple(alpha0 / sum_alpha, alpha1 / sum_alpha, alpha2 / sum_alpha);
    }
};

} // namespace PHARE::core

#endif
