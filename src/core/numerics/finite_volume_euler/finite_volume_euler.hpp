#ifndef PHARE_CORE_NUMERICS_FINITE_VOLUME_EULER_HPP
#define PHARE_CORE_NUMERICS_FINITE_VOLUME_EULER_HPP

#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/constants.hpp"
#include <iterator>
#include <tuple>

namespace PHARE::core
{
template<typename GridLayout>
class FiniteVolumeEuler_ref;

template<typename GridLayout>
class FiniteVolumeEuler : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename Field, typename... Fluxes>
    void operator()(Field const& U, Field& Unew, double const& dt, const Fluxes&... fluxes) const
    {
        if (!this->hasLayout())
            throw std::runtime_error("Error - FiniteVolumeEuler - GridLayout not set, cannot "
                                     "proceed to computation");

        FiniteVolumeEuler_ref{*this->layout_, dt}(U, Unew, fluxes...);
    }
};

template<typename GridLayout>
class FiniteVolumeEuler_ref
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    FiniteVolumeEuler_ref(GridLayout const& layout, double const dt)
        : layout_{layout}
        , dt_{dt}
    {
    }

    template<typename Field, typename... Fluxes>
    void operator()(Field const& U, Field& Unew, const Fluxes&... fluxes) const
    {
        layout_.evalOnBox(Unew, [&](auto&... args) mutable {
            finite_volume_euler_(U, Unew, {args...}, fluxes...);
        });
    }

private:
    GridLayout layout_;
    double dt_;

    template<typename Field, typename... Fluxes>
    void finite_volume_euler_(Field const& U, Field& Unew, MeshIndex<Field::dimension> index,
                              const Fluxes&... fluxes) const
    {
        auto&& flux_tuple = std::forward_as_tuple(fluxes...);

        auto&& F_x          = std::get<0>(flux_tuple);
        auto fluxCenteringX = layout_.centering(F_x.physicalQuantity());

        if constexpr (dimension == 1)
        {
            Unew(index)
                = U(index)
                  - (dt_ * layout_.inverseMeshSize(Direction::X))
                        * (F_x(layout_.nextIndex(fluxCenteringX[dirX], index[0])) - F_x(index));
        }
        else if constexpr (dimension >= 2)
        {
            auto&& F_y          = std::get<1>(flux_tuple);
            auto fluxCenteringY = layout_.centering(F_y.physicalQuantity());

            if constexpr (dimension == 2)
            {
                Unew(index)
                    = U(index)
                      - (dt_ * layout_.inverseMeshSize(Direction::X))
                            * (F_x(layout_.nextIndex(fluxCenteringX[dirX], index[0]), index[1])
                               - F_x(index))
                      - (dt_ * layout_.inverseMeshSize(Direction::Y))
                            * (F_y(index[0], layout_.nextIndex(fluxCenteringY[dirY], index[1]))
                               - F_y(index));
            }
            else if constexpr (dimension == 3)
            {
                auto&& F_z          = std::get<2>(flux_tuple);
                auto fluxCenteringZ = layout_.centering(F_z.physicalQuantity());

                Unew(index)
                    = U(index)
                      - (dt_ * layout_.inverseMeshSize(Direction::X))
                            * (F_x(layout_.nextIndex(fluxCenteringX[dirX], index[0]), index[1],
                                   index[2])
                               - F_x(index))
                      - (dt_ * layout_.inverseMeshSize(Direction::Y))
                            * (F_y(index[0], layout_.nextIndex(fluxCenteringY[dirY], index[1]),
                                   index[2])
                               - F_y(index))
                      - (dt_ * layout_.inverseMeshSize(Direction::Z))
                            * (F_z(index[0], index[1],
                                   layout_.nextIndex(fluxCenteringZ[dirZ], index[2]))
                               - F_z(index));
            }
        }
    }
};

} // namespace PHARE::core

#endif
