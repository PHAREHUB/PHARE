#ifndef PHARE_SRC_AMR_FIELD_FIELD_GEOMETRY_H
#define PHARE_SRC_AMR_FIELD_FIELD_GEOMETRY_H

#include <SAMRAI/hier/BoxGeometry.h>

#include "data/grid/gridlayoutdefs.h"
#include "field_overlap.h"
#include "utilities/types.h"

namespace PHARE
{
template<std::size_t dim, typename GridLayout, typename PhysicalQuantity>
class FieldGeometry : public SAMRAI::hier::BoxGeometry
{
public:
    /** \brief Construct a FieldGeometry on a region, for a specific quantity,
     * with a temporary gridlayout
     */
    FieldGeometry(SAMRAI::hier::Box const& box, GridLayout layout, PhysicalQuantity qty)
        : box_{box}
        , layout_{std::move(layout)}
        , quantity_{qty}
    {
    }




    /** \brief calculate overlap from two geometry boxes
     * When we want to calculate an overlap from two FieldGeometry, we give two boxes:
     * sourceMask : which is a region of the source patch and fillBox which is a region
     * on the destination patch, then we tell if we want to overwrite the interior of
     * the destination, after that we have to tell how the sourceMask should be transformed
     * to intersect it with the destination box (for example: when using periodics boundary,
     * the mask will be shifted).
     * Finaly we can also add some restrictions on the boxes for the destination.
     */
    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& destinationGeometry,
                     SAMRAI::hier::BoxGeometry const& sourceGeometry,
                     SAMRAI::hier::Box const& sourceMask, SAMRAI::hier::Box const& fillBox,
                     bool const overwriteInterior, SAMRAI::hier::Transformation const& sourceOffset,
                     bool const retry,
                     SAMRAI::hier::BoxContainer const& destinationRestrictBoxes
                     = SAMRAI::hier::BoxContainer{}) const final
    {
        std::shared_ptr<SAMRAI::hier::BoxOverlap> overlap;
        auto destinationCast = dynamic_cast<FieldGeometry const*>(&destinationGeometry);
        auto sourceCast      = dynamic_cast<FieldGeometry const*>(&sourceGeometry);

        if ((destinationCast != 0) && (sourceCast != 0))
        {
            overlap = doOverlap_(*destinationCast, *sourceCast, sourceMask, fillBox,
                                 overwriteInterior, sourceOffset, destinationRestrictBoxes);
        }
        else if (retry)
        {
            overlap = sourceGeometry.calculateOverlap(
                destinationGeometry, sourceGeometry, sourceMask, fillBox, overwriteInterior,
                sourceOffset, false, destinationRestrictBoxes);
        }
        return overlap;
    }




    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    setUpOverlap(SAMRAI::hier::BoxContainer const& boxes,
                 SAMRAI::hier::Transformation const& offset) const final
    {
        SAMRAI::hier::BoxContainer destinationBox;

        for (auto& box : boxes)
        {
            GridLayout layout = layoutFromBox_(box, layout_);
            SAMRAI::hier::Box fieldBox(toFieldBox(box, quantity_, layout, false));
            destinationBox.push_back(fieldBox);
        }

        return std::make_shared<FieldOverlap<dim>>(destinationBox, offset);
    }




    static SAMRAI::hier::Box toFieldBox(SAMRAI::hier::Box box, PhysicalQuantity qty,
                                        GridLayout const& layout, bool withGhost = true)
    {
        SAMRAI::hier::IntVector lower = box.lower();
        constexpr std::uint32_t dirX{0};
        constexpr std::uint32_t dirY{1};
        constexpr std::uint32_t dirZ{2};

        if (withGhost)
        {
            auto const& centering = GridLayout::centering(qty);

            SAMRAI::hier::IntVector shift(box.getDim());
            shift[dirX] = layout.nbrGhostNodes(centering[dirX]);

            if (dim > 1)
            {
                shift[dirY] = layout.nbrGhostNodes(centering[dirY]);
            }
            if (dim > 2)
            {
                shift[dirZ] = layout.nbrGhostNodes(centering[dirZ]);
            }

            lower = lower - shift;

            uint32 xStart = layout.ghostStartIndex(qty, Direction::X);
            uint32 xEnd   = layout.ghostEndIndex(qty, Direction::X);

            box.setLower(dirX, lower[dirX]);
            box.setUpper(dirX, xEnd - xStart + lower[dirX]);

            if (dim > 1)
            {
                uint32 yStart = layout.ghostStartIndex(qty, Direction::Y);
                uint32 yEnd   = layout.ghostEndIndex(qty, Direction::Y);

                box.setLower(dirY, lower[dirY]);
                box.setUpper(dirY, yEnd - yStart + lower[dirY]);
            }
            if (dim > 2)
            {
                uint32 zStart = layout.ghostStartIndex(qty, Direction::Z);
                uint32 zEnd   = layout.ghostEndIndex(qty, Direction::Z);

                box.setLower(dirZ, lower[dirZ]);
                box.setUpper(dirZ, zEnd - zStart + lower[dirZ]);
            }
        }
        else
        {
            uint32 xStart = layout.physicalStartIndex(qty, Direction::X);
            uint32 xEnd   = layout.physicalEndIndex(qty, Direction::X);

            box.setLower(dirX, lower[dirX]);
            box.setUpper(dirX, xEnd - xStart + lower[dirX]);

            if (dim > 1)
            {
                uint32 yStart = layout.physicalStartIndex(qty, Direction::Y);
                uint32 yEnd   = layout.physicalEndIndex(qty, Direction::Y);

                box.setLower(dirY, lower[dirY]);
                box.setUpper(dirY, yEnd - yStart + lower[dirY]);
            }
            if (dim > 2)
            {
                uint32 zStart = layout.physicalStartIndex(qty, Direction::Z);
                uint32 zEnd   = layout.physicalEndIndex(qty, Direction::Z);

                box.setLower(dirZ, lower[dirZ]);
                box.setUpper(dirZ, zEnd - zStart + lower[dirZ]);
            }
        }

        return box;
    }



private:
    SAMRAI::hier::Box box_;
    GridLayout layout_;
    PhysicalQuantity quantity_;



    void computeDestinationBoxes_(SAMRAI::hier::BoxContainer& destinationBoxes,
                                  FieldGeometry const& sourceGeometry,
                                  SAMRAI::hier::Box const& sourceMask,
                                  SAMRAI::hier::Box const& fillBox, bool const overwriteInterior,
                                  SAMRAI::hier::Transformation const& sourceOffset,
                                  SAMRAI::hier::BoxContainer const& destinationRestrictBoxes
                                  = SAMRAI::hier::BoxContainer()) const
    {
        SAMRAI::hier::Box sourceShift = sourceGeometry.box_ * sourceMask;

        GridLayout sourceShiftLayout = layoutFromBox_(sourceShift, sourceGeometry.layout_);
        GridLayout fillBoxLayout     = layoutFromBox_(fillBox, layout_);

        sourceOffset.transform(sourceShift);

        SAMRAI::hier::Box const destinationField(toFieldBox(box_, quantity_, layout_));


        // When we use an overlap, we want to copy data from the interior of source,
        // to the ghost of the destination, and if overwriteInterior is true , to the
        // interior;

        SAMRAI::hier::Box const sourceField(
            toFieldBox(sourceShift, quantity_, sourceShiftLayout, false));
        SAMRAI::hier::Box const fillField(toFieldBox(fillBox, quantity_, fillBoxLayout));
        SAMRAI::hier::Box const together(destinationField * sourceField * fillField);

        if (!together.empty())
        {
            if (!overwriteInterior)
            {
                SAMRAI::hier::Box interiorBox(toFieldBox(box_, quantity_, layout_, false));
                destinationBoxes.removeIntersections(together, interiorBox);
            }
            else
            {
                destinationBoxes.push_back(together);
            }
        }

        if (!destinationRestrictBoxes.empty() && !destinationBoxes.empty())
        {
            SAMRAI::hier::BoxContainer restrictBoxes;
            for (auto box = destinationRestrictBoxes.begin(); box != destinationRestrictBoxes.end();
                 ++box)
            {
                restrictBoxes.push_back(toFieldBox(*box, quantity_, layoutFromBox_(*box, layout_)));
            }
            destinationBoxes.intersectBoxes(restrictBoxes);
        }
    }




    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    doOverlap_(FieldGeometry const& destinationGeometry, FieldGeometry const& sourceGeometry,
               SAMRAI::hier::Box const& sourceMask, SAMRAI::hier::Box const& fillBox,
               bool const overwriteInterior, SAMRAI::hier::Transformation const& sourceOffset,
               SAMRAI::hier::BoxContainer const& destinationRestrictBoxes) const
    {
        SAMRAI::hier::BoxContainer destinationBox;

        destinationGeometry.computeDestinationBoxes_(destinationBox, sourceGeometry, sourceMask,
                                                     fillBox, overwriteInterior, sourceOffset,
                                                     destinationRestrictBoxes);

        return std::make_shared<FieldOverlap<dim>>(destinationBox, sourceOffset);
    }



    GridLayout layoutFromBox_(SAMRAI::hier::Box const& box, GridLayout const& layout) const
    {
        std::array<uint32, dim> nbCell;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            nbCell[iDim] = static_cast<uint32>(box.numberCells(iDim));
        }

        return GridLayout(layout.dxdydz(), nbCell, layout.layoutName(), layout.origin(),
                          layout.order());
    }
};




} // namespace PHARE

#endif
