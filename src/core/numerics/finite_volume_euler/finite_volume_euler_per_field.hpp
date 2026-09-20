#ifndef PHARE_CORE_NUMERICS_FINITE_VOLUME_EULER_HPP
#define PHARE_CORE_NUMERICS_FINITE_VOLUME_EULER_HPP


#include "core/utilities/constants.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include <tuple>


namespace PHARE::core
{


template<typename GridLayout>
class FiniteVolumeEulerPerField
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    FiniteVolumeEulerPerField(GridLayout const& layout, double const dt)
        : layout_{layout}
        , dt_{dt}
    {
    }

    template<typename Field, typename... Fluxes>
    void operator()(Field const& U, Field& Unew, Fluxes const&... fluxes) const
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
                              Fluxes const&... fluxes) const
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


// FiniteVolumeEulerPerField itself assumes plain (non-tiled) fields, since it needs a
// GridLayout matching the extent of the fields it is given (see evalOnBox). This wraps it so
// callers (see FiniteVolumeEuler) can pass either plain or tiled U/Unew/fluxes transparently.
template<typename GridLayout>
class FiniteVolumeEulerPerFieldSingleTransformer
{
public:
    FiniteVolumeEulerPerFieldSingleTransformer(GridLayout const& layout, double const dt)
        : layout_{layout}
        , dt_{dt}
    {
    }

    template<typename Field, typename... Fluxes>
    void operator()(Field const& U, Field& Unew, Fluxes const&... fluxes) const
    {
        if constexpr (is_field_tile_set_v<Field>)
        {
            auto& U_tiles    = U();
            auto& Unew_tiles = Unew();

            // U, Unew and each of fluxes are all defined on the same patch, so they share the
            // same tiling
            for (std::size_t i = 0; i < U_tiles.size(); ++i)
                FiniteVolumeEulerPerField<GridLayout>{U_tiles[i].layout(), dt_}(
                    U_tiles[i](), Unew_tiles[i](), fluxes()[i]()...);
            // probably need some inner sync here
        }
        else
            FiniteVolumeEulerPerField<GridLayout>{layout_, dt_}(U, Unew, fluxes...);
    }

private:
    GridLayout layout_;
    double dt_;
};

} // namespace PHARE::core

#endif
