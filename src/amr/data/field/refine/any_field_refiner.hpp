#ifndef PHARE_ANY_FIELD_REFINER_HPP
#define PHARE_ANY_FIELD_REFINER_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/data/grid/grid_tiles.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include "amr/utilities/box/amr_box.hpp"

#include <SAMRAI/hier/Box.h>

#include <cstddef>
#include <vector>
#include <optional>

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
            auto const coarseIdx = to_coarse(fineIndex);
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


    // same as operator() over every index of box, but for tiles each destination tile's
    // overlap and source tile candidates are resolved once, and one Derived is built per
    // tile pair, rather than scanning all tiles and building a Derived per fine index
    template<typename FieldT>
    void refine_box(FieldT const& coarseField, FieldT& fineField,
                    core::Box<int, dimension> const& box)
    {
        if constexpr (core::is_field_tile_set_v<FieldT>)
        {
            using SrcTile = std::decay_t<decltype(coarseField()[0])>;
            struct Candidate
            {
                SrcTile const* tile;
                core::Box<int, dimension> box; // field or ghost box, depending on the list
                std::optional<Derived> refiner{};
            };

            std::vector<Candidate> in_field, in_ghost;

            for (auto& dst_tile : fineField())
            {
                auto const dst_box = dst_tile.ghost_box();
                auto const overlap = box * dst_box;
                if (!overlap)
                    continue;

                // tile order is kept, so the first match per fine index is the same src tile
                // the full scan in operator() picks
                auto const coarse_overlap = core::Box<int, dimension>{to_coarse(overlap->lower),
                                                                      to_coarse(overlap->upper)};
                in_field.clear();
                in_ghost.clear();
                for (auto const& src_tile : coarseField())
                {
                    if (auto const field_box = src_tile.field_box(); field_box * coarse_overlap)
                        in_field.push_back({&src_tile, field_box});
                    if (auto const ghost_box = src_tile.ghost_box(); ghost_box * coarse_overlap)
                        in_ghost.push_back({&src_tile, ghost_box});
                }

                auto const refine_with = [&](Candidate& c, auto const& fineIndex) {
                    if (!c.refiner)
                        c.refiner.emplace(centerings_, samrai_box_from(dst_box),
                                          samrai_box_from(c.tile->ghost_box()), ratio_);
                    (*c.refiner)((*c.tile)(), dst_tile(), fineIndex);
                };
                auto const first_in = [](auto& candidates, auto const& coarseIdx) -> Candidate* {
                    for (auto& c : candidates)
                        if (isIn(coarseIdx, c.box))
                            return &c;
                    return nullptr;
                };

                for (auto const& fineIndex : *overlap)
                {
                    auto const coarseIdx = to_coarse(fineIndex);
                    if (auto* c = first_in(in_field, coarseIdx))
                        refine_with(*c, fineIndex);
                    else if (auto* c = first_in(in_ghost, coarseIdx))
                        refine_with(*c, fineIndex);
                }
            }
        }
        else
            for (auto const& fineIndex : box)
                static_cast<Derived&>(*this).refine(coarseField, fineField, fineIndex);
    }


protected:
    static auto to_coarse(core::Point<int, dimension> fineIndex)
    {
        for (auto& idx : fineIndex)
            idx = idx / refinementRatio;
        return fineIndex;
    }

    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
    std::array<core::QtyCentering, dimension> const centerings_;
    SAMRAI::hier::IntVector const& ratio_;
};
} // namespace PHARE::amr


#endif // PHARE_ANY_FIELD_REFINER_HPP
