#ifndef PHARE_AMR_UTILS_HPP
#define PHARE_AMR_UTILS_HPP

#include "core/def/phare_mpi.hpp"

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchData.h>


#include "amr/types/amr_types.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"
#include "core/def.hpp"

#include "amr/utilities/box/amr_box.hpp"

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
    void AMRToLocal(SAMRAI::hier::Box& AMRBox, SAMRAI::hier::Box const& referenceAMRBox);



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
        int constexpr dimension = GridLayoutT::dimension;

        SAMRAI::tbox::Dimension const dim{dimension};

        //  We get geometry information from the patch, such as meshSize, and physical origin
        auto patchGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            patch.getPatchGeometry());
        core::Point<double, dimension> origin;

        std::array<double, dimension> dl;

        if (patchGeom != nullptr)
        {
            auto pOrigin = patchGeom->getXLower();
            auto pDl     = patchGeom->getDx();

            for (std::size_t iDim = 0; iDim < dimension; ++iDim)
            {
                origin[iDim] = pOrigin[iDim];
                dl[iDim]     = pDl[iDim];
            }
        }
        else
        /*
          We assume that this is a temporary patch used by SAMRAI for data transfers
          Temporary patches are not given a Geometry at this moment so we can't use it.
          This happens in:
           SAMRAI::xfer::RefineTimeTransaction::packStream(tbox::MessageStream&stream)

          SEE: https://github.com/LLNL/SAMRAI/issues/147
        */
        {
            for (std::size_t iDim = 0; iDim < dimension; ++iDim)
            {
                origin[iDim] = 0;
                dl[iDim]     = 1;
            }
        }

        SAMRAI::hier::Box domain = patch.getBox();

        std::array<std::uint32_t, dimension> nbrCell;

        for (std::size_t iDim = 0; iDim < dimension; ++iDim)
        {
            nbrCell[iDim] = static_cast<std::uint32_t>(domain.numberCells(iDim));
        }

        return GridLayoutT{dl, nbrCell, origin, amr::Box<int, dimension>{domain}};
    }

    inline auto to_string(SAMRAI::hier::GlobalId const& id)
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
                                   std::forward<Action>(action), std::forward<Args...>(args...));
        }
    }

} // namespace amr
} // namespace PHARE

#endif // UTILS_HPP
