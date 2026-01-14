#ifndef CORE_NUMERICS_RECONSTRUCTION_MP5_HPP
#define CORE_NUMERICS_RECONSTRUCTION_MP5_HPP

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/slope_limiters/min_mod.hpp"
#include "core/utilities/index/index.hpp"
#include <utility>

namespace PHARE::core
{
template<typename GridLayout, typename SlopeLimiter = void>
class MP5Reconstruction
{
public:
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

        return std::make_pair(recons_mp5_L_(u_3, u_2, u_1, u, u1),
                              recons_mp5_R_(u_2, u_1, u, u1, u2));
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

        return std::make_pair(recons_mp5_L_(u_3, u_2, u_1, u, u1),
                              recons_mp5_R_(u_2, u_1, u, u1, u2));
    }

private:
    static auto recons_mp5_L_(auto ull, auto ul, auto u, auto ur, auto urr)
    {
        return recons_mp5_impl_(ull, ul, u, ur, urr);
    }

    static auto recons_mp5_R_(auto ull, auto ul, auto u, auto ur, auto urr)
    {
        return recons_mp5_impl_(urr, ur, u, ul, ull);
    }

    static auto recons_mp5_impl_(auto const v_m2, auto const v_m1, auto const u, auto const v_p1,
                                 auto const v_p2)
    {
        static constexpr auto alpha = 4.;
        auto const fi1_2  = (2. * v_m2 - 13. * v_m1 + 47. * u + 27. * v_p1 - 3. * v_p2) / 60.;
        auto const Dil    = u - v_m1;
        auto const Dir    = v_p1 - u;
        auto const fMP    = u + MinModLimiter::limit(Dir, alpha * Dil);
        auto const fUL    = u + alpha * Dil;
        auto const di     = Dir - Dil;
        auto const dir    = (v_p2 - v_p1) - Dir;
        auto const d4i1_2 = MinModLimiter::limit(4. * di - dir, 4. * dir - di, di, dir);
        auto const fMD    = (u + v_p1) / 2. - d4i1_2 / 2.;
        auto const fLC    = u + Dil / 2. + (4. / 3.) * d4i1_2;
        auto const fmin   = std::max(std::min({u, v_p1, fMD}), std::min({u, fUL, fLC}));
        auto const fmax   = std::min(std::max({u, v_p1, fMD}), std::max({u, fUL, fLC}));
        return (fi1_2 - u) * (fi1_2 - fMP) < 0.0 ? fi1_2 : std::clamp(fi1_2, fmin, fmax);
    }
};

} // namespace PHARE::core

#endif
