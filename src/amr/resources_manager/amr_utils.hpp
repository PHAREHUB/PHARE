#ifndef PHARE_AMR_UTILS_HPP
#define PHARE_AMR_UTILS_HPP


#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/constants.hpp"

#include "amr/types/amr_types.hpp"
#include "amr/utilities/box/amr_box.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/HierarchyNeighbors.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>

namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;

    /**
     * @brief offsetIsZero_ returns true of the transformation has zero offset
     */
    NO_DISCARD bool offsetIsZero(SAMRAI::hier::Transformation const& transformation);



    /**
     * @brief isSameBlock returns true if the transformation is not changing block
     */
    NO_DISCARD bool isSameBlock(SAMRAI::hier::Transformation const& transformation);



    /**
     * @brief AMRToLocal sets the AMRBox to local indexing relative to the referenceAMRBox
     */
    SAMRAI::hier::Box& AMRToLocal(SAMRAI::hier::Box& AMRBox,
                                  SAMRAI::hier::Box const& referenceAMRBox);



    /**
     * @brief AMRToLocal returns a local indexed box relative to the referenceAMRBox from the AMRBox
     */
    NO_DISCARD SAMRAI::hier::Box AMRToLocal(SAMRAI::hier::Box const& AMRBox,
                                            SAMRAI::hier::Box const& referenceAMRBox);



    /**
     * @brief AMRToLocal returns the vector to add to a box to put it in the local index space
     * relative to the referenceAMRBox
     */
    NO_DISCARD SAMRAI::hier::IntVector AMRToLocal(SAMRAI::hier::Box const& referenceAMRBox);


    /**
     * @brief localToAMR returns the vector to add to a box to put it in AMR index space from a
     * local index relative to referenceAMRBox
     */
    NO_DISCARD SAMRAI::hier::IntVector localToAMRVector(SAMRAI::hier::Box const& referenceAMRBox);


    /**
     * @brief AMRToLocal returns a local index relative to the referenceAMRBox lower bound
     *
     */
    template<std::size_t dimension, template<typename, std::size_t> typename Index>
    NO_DISCARD Index<int, dimension> AMRToLocal(Index<int, dimension> index,
                                                SAMRAI::hier::Box const& referenceAMRBox)
    {
        index[dirX] = index[dirX] - referenceAMRBox.lower(dirX);
        if constexpr (dimension > 1)
        {
            index[dirY] = index[dirY] - referenceAMRBox.lower(dirY);
        }
        if constexpr (dimension > 2)
        {
            index[dirZ] = index[dirZ] - referenceAMRBox.lower(dirZ);
        }
        return index;
    }
    /**
     * @brief localToAMR returns a amr index from a  relative index to the referenceAMRBox lower
     * bound
     *
     */
    template<std::size_t dimension, template<typename, std::size_t> typename Index>
    NO_DISCARD Index<int, dimension> localToAMR(Index<int, dimension> index,
                                                SAMRAI::hier::Box const& referenceAMRBox)
    {
        index[dirX] = index[dirX] + referenceAMRBox.lower(dirX);
        if constexpr (dimension > 1)
        {
            index[dirY] = index[dirY] + referenceAMRBox.lower(dirY);
        }
        if constexpr (dimension > 2)
        {
            index[dirZ] = index[dirZ] + referenceAMRBox.lower(dirZ);
        }
        return index;
    }

    /**
     * @brief toCoarseIndex returns a coarse index from a fine index assuming a ratio of 2 in AMR
     * index
     *
     * For positive indices, it is simply index / 2. For negative indices, index / 2 doesnt work
     * because of integer division. For example, -3 / 2 = -1, but we want it to be -2, while for
     * even negative integers e.g. -2 / 2 = -1, this is correct. We thus end up with the following
     * formula for negative indices: index / 2 + index % 2 (remember that in C++ the % is signed)
     * For example, -3 % 2 = -1, so we have -3 / 2 + -1 = -2.
     *
     */
    template<std::size_t dimension, template<typename, std::size_t> typename Index>
    NO_DISCARD Index<int, dimension> toCoarseIndex(Index<int, dimension> fineIndex)
    {
        auto coarseIdx{fineIndex};
        for (auto& idx : coarseIdx)
            idx = (idx >= 0) ? idx / 2 : idx / 2 + idx % 2;
        return coarseIdx;
    }

    /**
     * @brief refinedPosition returns an index refined index with the given ratio
     * bound
     *
     */
    template<std::size_t dimension, template<typename, std::size_t> typename Index>
    NO_DISCARD Index<int, dimension> refinedPosition(Index<int, dimension> index,
                                                     SAMRAI::hier::IntVector const& ratio)
    {
        index[dirX] *= ratio(dirX);
        if constexpr (dimension > 1)
        {
            index[dirY] *= ratio(dirY);
        }
        if constexpr (dimension > 2)
        {
            index[dirZ] *= ratio(dirZ);
        }
        return index;
    }


    template<typename GridLayoutT>
    NO_DISCARD GridLayoutT layoutFromPatch(SAMRAI::hier::Patch const& patch)
    {
        auto constexpr dimension = GridLayoutT::dimension;

        //  We get geometry information from the patch, such as meshSize, and physical origin
        auto patchGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            patch.getPatchGeometry());

        /*
          if(patchGeom == nullptr)

          We assume that this is a temporary patch used by SAMRAI for data transfers
          Temporary patches are not given a Geometry at this moment so we can't use it.
          This happens in:
           SAMRAI::xfer::RefineTimeTransaction::packStream(tbox::MessageStream&stream)

          SEE: https://github.com/LLNL/SAMRAI/issues/147
        */

        auto const dl = core::for_N_make_array<dimension>([&](auto i) {
            if (patchGeom != nullptr)
                return patchGeom->getDx()[i];
            return 1.0;
        });

        auto const origin = core::for_N_make_array<dimension>([&](auto i) {
            if (patchGeom != nullptr)
                return patchGeom->getXLower()[i];
            return 0.0;
        });

        auto const nbrCell = core::for_N_make_array<dimension>(
            [&](auto i) { return static_cast<std::uint32_t>(patch.getBox().numberCells(i)); });

        return {
            dl,
            nbrCell,
            origin, //
            amr::Box<int, dimension>{patch.getBox()},
            patch.getPatchLevelNumber(), //
        };
    }


    // potentially to replace with SAMRAI coarse to fine boundary stuff
    template<typename GridLayoutT> // fow now it gives us a box for only patch ghost layer
    NO_DISCARD auto makeNonLevelGhostBoxFor(SAMRAI::hier::Patch const& patch,
                                            SAMRAI::hier::PatchHierarchy const& hierarchy)
    {
        auto constexpr dimension       = GridLayoutT::dimension;
        auto const lvlNbr              = patch.getPatchLevelNumber();
        SAMRAI::hier::Box const domain = patch.getBox();
        auto const domBox              = phare_box_from<dimension>(domain);
        auto const particleGhostBox    = grow(domBox, GridLayoutT::nbrParticleGhosts());

        SAMRAI::hier::HierarchyNeighbors const hier_nbrs{hierarchy, lvlNbr, lvlNbr};
        auto const neighbors = hier_nbrs.getSameLevelNeighbors(domain, lvlNbr);
        std::vector<core::Box<int, GridLayoutT::dimension>> patchGhostLayerBoxes;
        patchGhostLayerBoxes.reserve(neighbors.size() + 1);
        patchGhostLayerBoxes.emplace_back(domBox);
        for (auto const& neighbox : neighbors)
            patchGhostLayerBoxes.emplace_back(
                *(particleGhostBox * phare_box_from<dimension>(neighbox)));

        return patchGhostLayerBoxes;
    }

    inline auto to_string(auto const& id)
    {
        std::stringstream patchID;
        patchID << id;
        return patchID.str();
    }

    template<typename GridLayout, typename ResMan, typename Action, typename... Args>
    void visitLevel(SAMRAI_Types::level_t& level, ResMan& resman, Action&& action, Args&&... args)
    {
        for (auto& patch : level)
        {
            auto guard        = resman.setOnPatch(*patch, args...);
            GridLayout layout = layoutFromPatch<GridLayout>(*patch);
            action(layout, to_string(patch->getGlobalId()),
                   static_cast<std::size_t>(level.getLevelNumber()));
        }
    }



    template<typename GridLayout, typename ResMan, typename Action, typename... Args>
    void visitHierarchy(SAMRAI::hier::PatchHierarchy& hierarchy, ResMan& resman, Action&& action,
                        int minLevel, int maxLevel, Args&&... args)
    {
        for (int iLevel = minLevel; iLevel < hierarchy.getNumberOfLevels() && iLevel <= maxLevel;
             iLevel++)
        {
            visitLevel<GridLayout>(*hierarchy.getPatchLevel(iLevel), resman,
                                   std::forward<Action>(action), std::forward<Args>(args)...);
        }
    }

} // namespace amr
} // namespace PHARE

#endif // UTILS_HPP
