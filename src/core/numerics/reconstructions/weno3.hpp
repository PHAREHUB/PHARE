#ifndef CORE_NUMERICS_RECONSTRUCTION_WENO3_HPP
#define CORE_NUMERICS_RECONSTRUCTION_WENO3_HPP

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include <utility>

namespace PHARE::core
{
template<typename GridLayout, typename SlopeLimiter = void>
class WENO3Reconstruction
{
public:
    using GridLayout_t = GridLayout;

    template<auto direction, typename Field>
    static auto reconstruct(Field const& F, MeshIndex<Field::dimension> index)
    {
        auto u_2 = F(GridLayout::template previous<direction>(
            GridLayout::template previous<direction>(index)));
        auto u_1 = F(GridLayout::template previous<direction>(index));
        auto u   = F(index);
        auto u1  = F(GridLayout::template next<direction>(index));

        return std::make_pair(recons_weno3_L_(u_2, u_1, u), recons_weno3_R_(u_1, u, u1));
    }

    template<auto direction, typename Field>
    static auto center_reconstruct(Field const& U, MeshIndex<Field::dimension> index,
                                   auto projection)
    {
        auto u_2 = GridLayout::project(U,
                                       GridLayout::template previous<direction>(
                                           GridLayout::template previous<direction>(index)),
                                       projection);
        auto u_1
            = GridLayout::project(U, GridLayout::template previous<direction>(index), projection);
        auto u  = GridLayout::project(U, index, projection);
        auto u1 = GridLayout::project(U, GridLayout::template next<direction>(index), projection);

        return std::make_pair(recons_weno3_L_(u_2, u_1, u), recons_weno3_R_(u_1, u, u1));
    }

private:
    static auto recons_weno3_L_(auto ul, auto u, auto ur)
    {
        static constexpr auto dL0 = 1. / 3.;
        static constexpr auto dL1 = 2. / 3.;

        auto const [wL0, wL1] = compute_weno3_weights(ul, u, ur, dL0, dL1);

        return wL0 * (-0.5 * ul + 1.5 * u) + wL1 * (0.5 * u + 0.5 * ur);
    }

    static auto recons_weno3_R_(auto ul, auto u, auto ur)
    {
        static constexpr auto dR0 = 2. / 3.;
        static constexpr auto dR1 = 1. / 3.;

        auto const [wR0, wR1] = compute_weno3_weights(ul, u, ur, dR0, dR1);

        return wR0 * (0.5 * u + 0.5 * ul) + wR1 * (-0.5 * ur + 1.5 * u);
    }

    static auto compute_weno3_weights(auto const ul, auto const u, auto const ur, auto const d0,
                                      auto const d1)
    {
        static constexpr auto eps = 1.e-6;

        auto const beta0 = (u - ul) * (u - ul);
        auto const beta1 = (ur - u) * (ur - u);

        auto const alpha0 = d0 / ((beta0 + eps) * (beta0 + eps));
        auto const alpha1 = d1 / ((beta1 + eps) * (beta1 + eps));

        auto const sum_alpha = alpha0 + alpha1;

        return std::make_tuple(alpha0 / sum_alpha, alpha1 / sum_alpha);
    }
};

} // namespace PHARE::core

#endif
