#ifndef CORE_NUMERICS_RECONSTRUCTION_LINEAR_HPP
#define CORE_NUMERICS_RECONSTRUCTION_LINEAR_HPP

#include "core/numerics/slope_limiters/van_leer.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include <utility>

namespace PHARE::core
{
template<typename GridLayout, typename SlopeLimiter = VanLeerLimiter>
class LinearReconstruction
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

        return std::make_pair(recons_linear_L_(u_2, u_1, u), recons_linear_R_(u_1, u, u1));
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

        return std::make_pair(recons_linear_L_(u_2, u_1, u), recons_linear_R_(u_1, u, u1));
    }

private:
    static auto recons_linear_L_(auto ul, auto u, auto ur)
    {
        auto Dil = (u - ul);
        auto Dir = (ur - u);
        auto Di  = SlopeLimiter::limit(Dil, Dir);

        return u + 0.5 * Di;
    }

    static auto recons_linear_R_(auto ul, auto u, auto ur)
    {
        auto Dil = (u - ul);
        auto Dir = (ur - u);
        auto Di  = SlopeLimiter::limit(Dil, Dir);

        return u - 0.5 * Di;
    }
};

} // namespace PHARE::core

#endif
