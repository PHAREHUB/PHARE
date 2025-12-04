#ifndef PHARE_SRC_AMR_FIELD_FIELD_GEOMETRY_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_GEOMETRY_HPP


#include "core/utilities/types.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include "field_overlap.hpp"

#include <SAMRAI/hier/Box.h>
#include "SAMRAI/hier/IntVector.h"
#include <SAMRAI/hier/BoxGeometry.h>

#include <cassert>

namespace PHARE
{
namespace amr
{


    // this class is needed because the FieldFillPattern needs to cast the
    // generic BoxGeometry into the specific geometry but cannot cast into
    // the FieldGeometry below because it does not have the GridLayoutT and
    // PhysicalQuantity for template arguments.
    template<std::size_t dimension>
    class FieldGeometryBase : public SAMRAI::hier::BoxGeometry
    {
    public:
        virtual ~FieldGeometryBase() {}
        FieldGeometryBase(SAMRAI::hier::Box const& patch_box,
                          SAMRAI::hier::Box const& ghostFieldBox,
                          SAMRAI::hier::Box const& interiorFieldBox,
                          std::array<core::QtyCentering, dimension> const& centerings)
            : patchBox{patch_box}
            , ghostFieldBox_{ghostFieldBox}
            , interiorFieldBox_{interiorFieldBox}
            , centerings_{centerings}
        {
        }

        auto const& interiorFieldBox() const { return interiorFieldBox_; }

        SAMRAI::hier::Box const patchBox;

    protected:
        SAMRAI::hier::Box const ghostFieldBox_;
        SAMRAI::hier::Box const interiorFieldBox_;
        std::array<core::QtyCentering, dimension> const centerings_;
    };

    template<typename GridLayoutT, typename PhysicalQuantity>
    /**
     * @brief The FieldGeometry class
     */
    class FieldGeometry : public FieldGeometryBase<GridLayoutT::dimension>
    {
    public:
        using Super                               = FieldGeometryBase<GridLayoutT::dimension>;
        static constexpr std::size_t dimension    = GridLayoutT::dimension;
        static constexpr std::size_t interp_order = GridLayoutT::interp_order;

        /** \brief Construct a FieldGeometry on a region, for a specific quantity,
         * with a temporary gridlayout
         */
        FieldGeometry(SAMRAI::hier::Box const& box, GridLayoutT const& layout,
                      PhysicalQuantity const qty)
            : Super(box,
                    toFieldBox(SAMRAI::hier::Box::grow(
                                   box, SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dimension},
                                                                GridLayoutT::nbrGhosts()}),
                               qty, layout),
                    toFieldBox(box, qty, layout), GridLayoutT::centering(qty))
            , layout_{layout}
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
         * Finally we can also add some restrictions on the boxes for the destination.
         */
        std::shared_ptr<SAMRAI::hier::BoxOverlap> calculateOverlap(
            SAMRAI::hier::BoxGeometry const& destinationGeometry,
            SAMRAI::hier::BoxGeometry const& sourceGeometry, SAMRAI::hier::Box const& sourceMask,
            SAMRAI::hier::Box const& fillBox, bool const overwriteInterior,
            SAMRAI::hier::Transformation const& sourceOffset, [[maybe_unused]] bool const retry,
            SAMRAI::hier::BoxContainer const& destinationRestrictBoxes
            = SAMRAI::hier::BoxContainer{}) const final
        {
            auto& destinationCast = dynamic_cast<FieldGeometry const&>(destinationGeometry);
            auto& sourceCast      = dynamic_cast<FieldGeometry const&>(sourceGeometry);
            return doOverlap_(destinationCast, sourceCast, sourceMask, fillBox, overwriteInterior,
                              sourceOffset, destinationRestrictBoxes);
        }




        std::shared_ptr<SAMRAI::hier::BoxOverlap>
        setUpOverlap(SAMRAI::hier::BoxContainer const& boxes,
                     SAMRAI::hier::Transformation const& offset) const final
        {
            SAMRAI::hier::BoxContainer destinationBoxes;

            for (auto& box : boxes)
            {
                core::GridLayout const layout = layoutFromBox(box, layout_);
                destinationBoxes.push_back(toFieldBox(box, quantity_, layout));
            }

            return std::make_shared<FieldOverlap>(destinationBoxes, offset);
        }




        /**
         * @brief toFieldBox takes an AMR cell-centered box and creates a box
         * that is adequate for the specified quantity. The layout is used to know
         * the centering, nbr of ghosts of the specified quantity.
         * @param withGhost true if we want to include the ghost nodes in the field box.
         *
         * Note : precondition : the nbr of physical cells of the layout must correspond
         *  to the nbr of cells of the box.
         */
        static SAMRAI::hier::Box toFieldBox(SAMRAI::hier::Box box, PhysicalQuantity qty,
                                            GridLayoutT const& layout)
        {
            auto const lower = box.lower();
            auto const upper = box.upper();
            using PHARE::core::dirX;
            using PHARE::core::dirY;
            using PHARE::core::dirZ;
            // lower/upper of 'box' are AMR cell-centered coordinates

            // we want a box that starts at lower
            // the upper index is

            // example:
            // . = primal node
            // x = dual node
            // here we have nbrGhost = 1 for both

            //     0     1     2     3     4     5     6        (local dual Index)
            //  0     1     2     3     4     5     6     7     (local primal Index)
            //  .  x  |  x  .  x  .  x  .  x  .  x  |  x  .
            //           ^
            //         box.lower (AMR index)

            // for primal = AMR upper index will be box.lower + 6 (pei) -1(psi) = 5 without ghosts
            // for dual = AMR upper index will be  box.lower + 5(pei) - 1(psi) = lower+ 4 without
            // ghosts

            // with ghosts :
            // box.lower must be shifted left to move to the first ghost node
            // box.upper is still box.lower + end-start, end &start of ghosts

            auto const centerings = layout.centering(qty);
            core::for_N<dimension>( //
                [&](auto i) {
                    box.setLower(i, lower[i]);
                    auto const is_primal = (centerings[i] == core::QtyCentering::primal) ? 1 : 0;
                    box.setUpper(i, upper[i] + is_primal);
                } //
            );


            return box;
        }

        static SAMRAI::hier::BoxContainer toFieldBoxes(SAMRAI::hier::BoxContainer const& boxes,
                                                       PhysicalQuantity qty,
                                                       GridLayoutT const& layout)
        {
            auto to_field_boxes = boxes;
            for (auto& box : to_field_boxes)
                box = toFieldBox(box, qty, layout);
            return to_field_boxes;
        }

        /**
         * @brief The origin of the returned layout should NOT be used
         * this is only to get start and end index for physical and ghost
         *
         * Note : this function is used whenever we don't have a layout with the same
         * nbr of cells as the box (ex: we want restriction of a patch box etc.)
         */
        static GridLayoutT layoutFromBox(SAMRAI::hier::Box const& box, GridLayoutT const& layout)
        {
            std::array<std::uint32_t, dimension> nbCell;
            for (std::size_t iDim = 0; iDim < dimension; ++iDim)
            {
                nbCell[iDim] = static_cast<std::uint32_t>(box.numberCells(iDim));
            }

            return GridLayoutT(layout.meshSize(), nbCell, layout.origin());
        }


    private:
        GridLayoutT layout_;
        PhysicalQuantity quantity_;


        /*** \brief Compute destination box representing the intersection of two geometry
         *
         *   \param destinationBoxes BoxContainer that will be filled of box
         *   \param sourceGeometry represent the geometry of the source data
         *   \param sourceMask restrict the portion concerned by the source data
         *   \param fillBox restrict the portion where data will be put on the destination
         *   \param destinationRestrictBoxes container of box that will restrict the intersection
         *
         */
        void computeDestinationBoxes_(SAMRAI::hier::BoxContainer& destinationBoxes,
                                      FieldGeometry const& sourceGeometry,
                                      SAMRAI::hier::Box const& sourceMask,
                                      SAMRAI::hier::Box const& fillBox,
                                      bool const overwriteInterior,
                                      SAMRAI::hier::Transformation const& sourceOffset,
                                      SAMRAI::hier::BoxContainer const& destinationRestrictBoxes
                                      = SAMRAI::hier::BoxContainer()) const
        {
            // we have three boxes :
            // - the sourceBox : the is where the data is to be taken from
            // - the destinationBox : this is the whole box on the "destination"
            // - the fillBox: this is a restriction of the destinationBox

            // all these boxes are in AMR (cell-centered) indexes.
            // we need to translate them into boxes describing the proper field layout.
            // (for instance, a cell-centered box is not adequatly describing a node centered
            // quantity).

            // the destinationsBoxes is to be filled with:
            // the intersection of the sourceBox, the destinationBox and the fillBox
            // with potential restrictionBoxes.
            // In case we don't want to fill interior nodes, we remove the interior
            // of the destination box, from the above intersection, which may result
            // in adding multiple boxes into the container.


            // the sourceMask is a restriction of the sourceBox
            // so we need to intersect it with the sourceBox, then to apply a transformation
            // to account for the periodicity
            SAMRAI::hier::Box sourceShift = sourceGeometry.ghostFieldBox_ * sourceMask;
            sourceOffset.transform(sourceShift);


            // ok let's get the boxes for the fields from cell-centered boxes now
            bool withGhosts = true;

            core::GridLayout sourceShiftLayout = layoutFromBox(sourceShift, sourceGeometry.layout_);
            core::GridLayout fillBoxLayout     = layoutFromBox(fillBox, layout_);

            auto const& destinationBox = this->ghostFieldBox_;

            SAMRAI::hier::Box const sourceBox{
                toFieldBox(sourceShift, quantity_, sourceShiftLayout)};

            SAMRAI::hier::Box const fillField{toFieldBox(fillBox, quantity_, fillBoxLayout)};


            // now we have all boxes shifted and translated to field boxes
            // let's compute the interesection of them all.

            SAMRAI::hier::Box const together(destinationBox * sourceBox * fillField);

            // if the interesection is not empty we either push it into the container
            // if we don't want to fill the interior we remove it from the intersection
            // which may add multiple boxes to the container.

            if (!together.empty())
            {
                if (overwriteInterior)
                {
                    destinationBoxes.push_back(together);
                }
                else
                {
                    destinationBoxes.removeIntersections(together, this->interiorFieldBox_);
                }
            }

            if (!destinationRestrictBoxes.empty() && !destinationBoxes.empty())
            {
                SAMRAI::hier::BoxContainer restrictBoxes;
                for (auto box = destinationRestrictBoxes.begin();
                     box != destinationRestrictBoxes.end(); ++box)
                {
                    restrictBoxes.push_back(
                        toFieldBox(*box, quantity_, layoutFromBox(*box, layout_)));
                }

                // will only keep of together the boxes that interesect the restrictions
                destinationBoxes.intersectBoxes(restrictBoxes);
            }
        }



        /**
         * @brief doOverlap_ will return a field overlap from the source to the dest
         * geometry. This function makes this overlap by calculating the (possibly multiple)
         * destination box(es) and the transformation from source to dest.
         */
        std::shared_ptr<SAMRAI::hier::BoxOverlap>
        doOverlap_(FieldGeometry const& destinationGeometry, FieldGeometry const& sourceGeometry,
                   SAMRAI::hier::Box const& sourceMask, SAMRAI::hier::Box const& fillBox,
                   bool const overwriteInterior, SAMRAI::hier::Transformation const& sourceOffset,
                   SAMRAI::hier::BoxContainer const& destinationRestrictBoxes) const
        {
            SAMRAI::hier::BoxContainer destinationBoxes;

            destinationGeometry.computeDestinationBoxes_(destinationBoxes, sourceGeometry,
                                                         sourceMask, fillBox, overwriteInterior,
                                                         sourceOffset, destinationRestrictBoxes);

            return std::make_shared<FieldOverlap>(destinationBoxes, sourceOffset);
        }
    };

} // namespace amr


} // namespace PHARE

#endif
