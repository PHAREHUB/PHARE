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
        double eps = 1.e-6;

        double dL0 = 1. / 3.;
        double dL1 = 2. / 3.;

        auto beta0 = (u - ul) * (u - ul);
        auto beta1 = (ur - u) * (ur - u);

        auto alphaL0 = dL0 / ((beta0 + eps) * (beta0 + eps));
        auto alphaL1 = dL1 / ((beta1 + eps) * (beta1 + eps));

        auto wL0 = alphaL0 / (alphaL0 + alphaL1);
        auto wL1 = alphaL1 / (alphaL0 + alphaL1);

        return wL0 * (-0.5 * ul + 1.5 * u) + wL1 * (0.5 * u + 0.5 * ur);
    }

    static auto recons_weno3_R_(auto ul, auto u, auto ur)
    {
        double eps = 1.e-6;

        double dR0 = 2. / 3.;
        double dR1 = 1. / 3.;

        auto beta1 = (u - ul) * (u - ul);
        auto beta2 = (ur - u) * (ur - u);

        auto alphaR0 = dR0 / ((beta1 + eps) * (beta1 + eps));
        auto alphaR1 = dR1 / ((beta2 + eps) * (beta2 + eps));

        auto wR0 = alphaR0 / (alphaR0 + alphaR1);
        auto wR1 = alphaR1 / (alphaR0 + alphaR1);

        return wR0 * (0.5 * u + 0.5 * ul) + wR1 * (-0.5 * ur + 1.5 * u);
    }
};

} // namespace PHARE::core

#endif
