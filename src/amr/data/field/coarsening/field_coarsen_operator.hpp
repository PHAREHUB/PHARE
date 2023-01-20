#ifndef PHARE_FIELD_DATA_COARSEN_HPP
#define PHARE_FIELD_DATA_COARSEN_HPP

#include "amr/data/field/field_data.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "field_coarsen_index_weight.hpp"
#include "default_field_coarsener.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/IntVector.h>


namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    //
    template<typename GridLayoutT, typename FieldT, typename FieldCoarsenerPolicy,
             typename PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity())>
    /**
     * @brief The FieldCoarsenOperator class
     */
    class FieldCoarsenOperator : public SAMRAI::hier::CoarsenOperator
    {
        static constexpr std::size_t n_ghosts
            = GridLayoutT::template nbrGhosts<core::QtyCentering, core::QtyCentering::dual>();

    public:
        static constexpr std::size_t dimension = GridLayoutT::dimension;
        using FieldDataT                       = FieldData<GridLayoutT, FieldT>;

        FieldCoarsenOperator()
            : SAMRAI::hier::CoarsenOperator("FieldDataCoarsenOperator")
        {
        }

        FieldCoarsenOperator(FieldCoarsenOperator const&)            = delete;
        FieldCoarsenOperator(FieldCoarsenOperator&&)                 = delete;
        FieldCoarsenOperator& operator=(FieldCoarsenOperator const&) = delete;
        FieldCoarsenOperator&& operator=(FieldCoarsenOperator&&)     = delete;


        virtual ~FieldCoarsenOperator() = default;




        /** @brief return the priority of the operator
         *  this return 0, meaning that this operator
         * have the most priority
         */
        int getOperatorPriority() const override { return 0; }




        /** @brief Return the stencil width associated with the coarsening operator.
         *
         *  The SAMRAI transfer routines guarantee that the source patch will contain
         * sufficient ghostCell data surrounding the interior to satisfy the stencil
         * width requirements for each coarsening operator.
         *
         * In our case, we allow a RF up to 10, so having 5 ghost width is sufficient
         *
         */
        SAMRAI::hier::IntVector getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
        {
            return SAMRAI::hier::IntVector{dim, 2};
        }




        /** @brief given a coarseBox, coarse data from the fine patch on the intersection of this
         * box and the box of the destination (the box of the coarse patch).
         *
         * This method will extract fieldData from the two patches, and then
         * get the Field and GridLayout encapsulated into the fieldData.
         * With the help of FieldGeometry, transform the coarseBox to the correct index.
         * After that we can now create FieldCoarsen with the indexAndWeight implementation
         * selected. Finnaly loop over the indexes in the box, and apply the coarsening defined in
         * FieldCoarsen operator
         *
         */
        void coarsen(SAMRAI::hier::Patch& destinationPatch, SAMRAI::hier::Patch const& sourcePatch,
                     int const destinationId, int const sourceId,
                     SAMRAI::hier::Box const& coarseBox,
                     SAMRAI::hier::IntVector const& ratio) const override
        {
            auto& destinationField        = FieldDataT::getField(destinationPatch, destinationId);
            auto const& sourceField       = FieldDataT::getField(sourcePatch, sourceId);
            auto const& sourceLayout      = FieldDataT::getLayout(sourcePatch, sourceId);
            auto const& destinationLayout = FieldDataT::getLayout(destinationPatch, destinationId);

            // we assume that quantity are the same
            // note that an assertion will be raised
            // in coarseIt operator
            auto const& qty = destinationField.physicalQuantity();


            bool const withGhost{true};

            // We get different boxes : destination , source, restrictBoxes
            // and transform them in the correct indexing.
            auto destinationBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
                destinationPatch.getBox(), qty, destinationLayout, withGhost);

            auto sourceBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
                sourcePatch.getBox(), qty, sourceLayout, withGhost);

            auto coarseLayout = FieldGeometry<GridLayoutT, PhysicalQuantity>::layoutFromBox(
                coarseBox, destinationLayout);

            auto coarseFieldBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
                coarseBox, qty, coarseLayout, !withGhost);

            // finnaly we compute the intersection
            auto intersectionBox = destinationBox * coarseFieldBox;




            // We can now create the coarsening operator
            FieldCoarsenerPolicy coarsener{destinationLayout.centering(qty), sourceBox,
                                           destinationBox, ratio};

            // now we can loop over the intersection box

            core::Point<int, dimension> startIndex;
            core::Point<int, dimension> endIndex;

            startIndex[dirX] = intersectionBox.lower(dirX);
            endIndex[dirX]   = intersectionBox.upper(dirX);

            if constexpr (dimension > 1)
            {
                startIndex[dirY] = intersectionBox.lower(dirY);
                endIndex[dirY]   = intersectionBox.upper(dirY);
            }
            if constexpr (dimension > 2)
            {
                startIndex[dirZ] = intersectionBox.lower(dirZ);
                endIndex[dirZ]   = intersectionBox.upper(dirZ);
            }

            if constexpr (dimension == 1)
            {
                for (int ix = startIndex[dirX]; ix <= endIndex[dirX]; ++ix)
                {
                    coarsener(sourceField, destinationField, {{ix}});
                }
            }




            else if constexpr (dimension == 2)
            {
                for (int ix = startIndex[dirX]; ix <= endIndex[dirX]; ++ix)
                {
                    for (int iy = startIndex[dirY]; iy <= endIndex[dirY]; ++iy)
                    {
                        coarsener(sourceField, destinationField, {{ix, iy}});
                    }
                }
            }




            else if constexpr (dimension == 3)
            {
                for (int ix = startIndex[dirX]; ix <= endIndex[dirX]; ++ix)
                {
                    for (int iy = startIndex[dirY]; iy <= endIndex[dirY]; ++iy)
                    {
                        for (int iz = startIndex[dirZ]; iz <= endIndex[dirZ]; ++iz)

                        {
                            coarsener(sourceField, destinationField, {{ix, iy, iz}});
                        }
                    }
                }
            } // end 3D
        }
    };
} // namespace amr
} // namespace PHARE


#endif
