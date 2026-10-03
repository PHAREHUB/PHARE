#ifndef PHARE_ANY_FIELD_REFINER_HPP
#define PHARE_ANY_FIELD_REFINER_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/data/grid/grid_tiles.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include "amr/utilities/box/amr_box.hpp"

#include <SAMRAI/hier/Box.h>

#include <cstddef>

namespace PHARE::amr
{

/** \brief CRTP base factoring out the tiled-vs-plain dispatch shared by all point-wise field
 * refiners (ElectricFieldRefiner, MagneticFieldRefiner, MagneticFieldInitRefiner,
 * MagneticFieldRegrider, MHDFieldRefiner, MHDFluxRefiner, ...).
 *
 * Derived classes inherit the constructor and operator()(), and only need to implement
 * refine(FieldT const&, FieldT&, Point) for the plain (non tiled) case. When FieldT is a
 * FieldTileSet, operator() finds the owning source/destination tiles for fineIndex and
 * recurses on a Derived built from those tiles' own (smaller) ghost boxes, until FieldT
 * is plain and refine() is called.
 */
template<typename Derived, std::size_t dimension>
class AnyFieldRefiner
{
public:
    AnyFieldRefiner(auto const& centering, auto const& destinationGhostBox,
                    auto const& sourceGhostBox, auto const& ratio)
        : fineBox_{destinationGhostBox}
        , coarseBox_{sourceGhostBox}
        , centerings_{centering}
        , ratio_{ratio}
    {
    }


    template<typename FieldT>
    void operator()(FieldT const& coarseField, FieldT& fineField,
                    core::Point<int, dimension> fineIndex)
    {
        if constexpr (core::is_field_tile_set_v<FieldT>)
        {
            auto coarseIdx{fineIndex};
            for (auto& idx : coarseIdx)
                idx = idx / refinementRatio;
            for (auto& dst_tile : fineField())
                if (auto const dst_box = dst_tile.ghost_box(); isIn(fineIndex, dst_box))
                {
                    auto const do_refine = [&](auto const& src_tile) {
                        Derived{centerings_, samrai_box_from(dst_box),
                                samrai_box_from(src_tile.ghost_box()),
                                ratio_}(src_tile(), dst_tile(), fineIndex);
                    };
                    bool found = false;
                    for (auto const& src_tile : coarseField())
                        if (isIn(coarseIdx, src_tile.field_box()))
                        {
                            do_refine(src_tile);
                            found = true;
                            break;
                        }
                    if (!found)
                        for (auto const& src_tile : coarseField())
                            if (isIn(coarseIdx, src_tile.ghost_box()))
                            {
                                do_refine(src_tile);
                                break;
                            }
                }
        }
        else
            static_cast<Derived&>(*this).refine(coarseField, fineField, fineIndex);
    }


protected:
    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
    std::array<core::QtyCentering, dimension> const centerings_;
    SAMRAI::hier::IntVector const& ratio_;
};
} // namespace PHARE::amr


#endif // PHARE_ANY_FIELD_REFINER_HPP
