#ifndef PHARE_AMR_UTILS_H
#define PHARE_AMR_UTILS_H

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchData.h>

#include "types/amr_types.h"
#include "utilities/box/box.h"
#include "utilities/constants.h"
#include "utilities/point/point.h"




namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    using core::uint32;
    /**
     * @brief offsetIsZero_ returns true of the transformation has zero offset
     */
    bool offsetIsZero(SAMRAI::hier::Transformation const& transformation);



    /**
     * @brief isSameBlock returns true if the transformation is not changing block
     */
    bool isSameBlock(SAMRAI::hier::Transformation const& transformation);



    /**
     * @brief AMRToLocal sets the AMRBox to local indexing relative to the referenceAMRBox
     */
    void AMRToLocal(SAMRAI::hier::Box& AMRBox, SAMRAI::hier::Box const& referenceAMRBox);



    /**
     * @brief AMRToLocal returns a local indexed box relative to the referenceAMRBox from the AMRBox
     */
    SAMRAI::hier::Box AMRToLocal(SAMRAI::hier::Box const& AMRBox,
                                 SAMRAI::hier::Box const& referenceAMRBox);



    /**
     * @brief AMRToLocal returns the vector to add to a box to put it in the local index space
     * relative to the referenceAMRBox
     */
    SAMRAI::hier::IntVector AMRToLocal(SAMRAI::hier::Box const& referenceAMRBox);


    /**
     * @brief localToAMR returns the vector to add to a box to put it in AMR index space from a
     * local index relative to referenceAMRBox
     */
    SAMRAI::hier::IntVector localToAMRVector(SAMRAI::hier::Box const& referenceAMRBox);


    /**
     * @brief AMRToLocal returns a local index relative to the referenceAMRBox lower bound
     *
     */
    template<std::size_t dimension, template<typename, std::size_t> typename Index>
    Index<int, dimension> AMRToLocal(Index<int, dimension> index,
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
    Index<int, dimension> localToAMR(Index<int, dimension> index,
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
    Index<int, dimension> refinedPosition(Index<int, dimension> index,
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


    template<std::size_t dim>
    auto toPHAREBox(SAMRAI::hier::Box b)
    {
        if constexpr (dim == 1)
        {
            auto lower = core::Point{b.lower()[0]};
            auto upper = core::Point{b.upper()[0]};
            return core::Box{lower, upper};
        }

        else if constexpr (dim == 2)
        {
            auto lower = core::Point{b.lower()[0], b.lower()[1]};
            auto upper = core::Point{b.upper()[0], b.upper()[1]};
            return core::Box{lower, upper};
        }

        else if constexpr (dim == 3)
        {
            auto lower = core::Point{b.lower()[0], b.lower()[1], b.lower()[2]};
            auto upper = core::Point{b.upper()[0], b.upper()[1], b.upper()[2]};
            return core::Box{lower, upper};
        }
    }



    template<typename GridLayoutT>
    GridLayoutT layoutFromPatch(SAMRAI::hier::Patch const& patch)
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
        {
            // in case that the patch does not have a CartesianPatchGeometry
            // the gridlayout will most likely throw at the construction
            // so we may throw here instead
            throw std::runtime_error(
                "The geometry on the patch is not set, please verify your configuration");
        }

        SAMRAI::hier::Box domain = patch.getBox();

        std::array<uint32, dimension> nbrCell;

        for (std::size_t iDim = 0; iDim < dimension; ++iDim)
        {
            nbrCell[iDim] = static_cast<uint32>(domain.numberCells(iDim));
        }

        return GridLayoutT{dl, nbrCell, origin, toPHAREBox<dimension>(domain)};
    }




    template<typename GridLayout, typename ResMan, typename Action, typename... Args>
    void visitLevel(SAMRAI_Types::level_t& level, ResMan& resman, Action&& action, Args&&... args)
    {
        for (auto& patch : level)
        {
            auto guard        = resman.setOnPatch(*patch, args...);
            GridLayout layout = layoutFromPatch<GridLayout>(*patch);
            std::stringstream patchID;
            patchID << patch->getGlobalId();
            action(layout, patchID.str(), level.getLevelNumber());
        }
    }


    template<typename GridLayout, typename ResMan, typename Action, typename... Args>
    void visitHierarchy(SAMRAI::hier::PatchHierarchy& hierarchy, ResMan& resman, Action&& action,
                        int minLevel, int maxLevel, Args&&... args)
    {
        assert(hierarchy.getNumberOfLevels() > maxLevel);
        for (int iLevel = minLevel; iLevel <= maxLevel; iLevel++)
        {
            visitLevel<GridLayout>(*hierarchy.getPatchLevel(iLevel), resman, action, args...);
        }
    }

} // namespace amr
} // namespace PHARE

#endif // UTILS_H
