#ifndef PHARE_FIELD_DATA_COARSEN_HPP
#define PHARE_FIELD_DATA_COARSEN_HPP

#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"

#include "amr/data/field/field_data.hpp"
#include "amr/data/field/field_geometry.hpp"

#include "default_field_coarsener.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/IntVector.h>


namespace PHARE::amr
{


template<typename Dst>
void coarsen_field(Dst& destinationField, auto& sourceField, auto& intersectionBox, auto& coarsener)
{
    auto constexpr static dimension = Dst::dimension;

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

} // namespace PHARE::amr


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




        /** @brief given a coarseBox, coarse data from the fine patch on the intersection of
         * this box and the box of the destination (the box of the coarse patch).
         *
         * This method will extract fieldData from the two patches, and then
         * get the Field and GridLayout encapsulated into the fieldData.
         * With the help of FieldGeometry, transform the coarseBox to the correct index.
         * After that we can now create FieldCoarsen with the indexAndWeight implementation
         * selected. Finaly loop over the indexes in the box, and apply the coarsening defined
         * in FieldCoarsen operator
         *
         */
        void coarsen(SAMRAI::hier::Patch& destinationPatch, SAMRAI::hier::Patch const& sourcePatch,
                     int const destinationId, int const sourceId,
                     SAMRAI::hier::Box const& coarseBox,
                     SAMRAI::hier::IntVector const& ratio) const override
        {
            auto& destinationField   = FieldDataT::getField(destinationPatch, destinationId);
            auto const& sourceField  = FieldDataT::getField(sourcePatch, sourceId);
            auto const& sourceLayout = FieldDataT::getLayout(sourcePatch, sourceId);
            auto const& destLayout   = FieldDataT::getLayout(destinationPatch, destinationId);
            using FieldGeometryT     = FieldGeometry<GridLayoutT, PhysicalQuantity>;

            // we assume that quantity are the same
            // note that an assertion will be raised
            // in coarseIt operator
            auto const& qty = destinationField.physicalQuantity();

            // We get different boxes : destination , source, restrictBoxes
            // and transform them in the correct indexing.
            auto destPData = destinationPatch.getPatchData(destinationId);
            auto srcPData  = sourcePatch.getPatchData(sourceId);
            auto destGBox  = FieldGeometryT::toFieldBox(destPData->getGhostBox(), qty, destLayout);
            auto srcGBox   = FieldGeometryT::toFieldBox(srcPData->getGhostBox(), qty, sourceLayout);
            auto coarseLayout          = FieldGeometryT::layoutFromBox(coarseBox, destLayout);
            auto coarseFieldBox        = FieldGeometryT::toFieldBox(coarseBox, qty, coarseLayout);
            auto const intersectionBox = destGBox * coarseFieldBox;
            // We can now create the coarsening operator
            FieldCoarsenerPolicy coarsener{destLayout.centering(qty), srcGBox, destGBox, ratio};

            coarsen_field(destinationField, sourceField, intersectionBox, coarsener);
        }
    };
} // namespace amr
} // namespace PHARE


namespace PHARE::amr
{


template<std::size_t rank, typename GridLayoutT, typename FieldT, typename FieldCoarsenerPolicy,
         typename PhysicalQuantity>
class TensorFieldCoarsenOperator : public SAMRAI::hier::CoarsenOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using TensorFieldDataT = TensorFieldData<rank, GridLayoutT, FieldT, PhysicalQuantity>;
    using FieldDataT       = FieldData<GridLayoutT, FieldT>;

    static constexpr std::size_t N = TensorFieldDataT::N;

    TensorFieldCoarsenOperator()
        : SAMRAI::hier::CoarsenOperator("FieldDataCoarsenOperator")
    {
    }

    TensorFieldCoarsenOperator(TensorFieldCoarsenOperator const&)            = delete;
    TensorFieldCoarsenOperator(TensorFieldCoarsenOperator&&)                 = delete;
    TensorFieldCoarsenOperator& operator=(TensorFieldCoarsenOperator const&) = delete;
    TensorFieldCoarsenOperator&& operator=(TensorFieldCoarsenOperator&&)     = delete;


    virtual ~TensorFieldCoarsenOperator() = default;




    /** @brief return the priority of the operator
     *  this return 0, meaning that this operator have the most priority
     */
    int getOperatorPriority() const override { return 0; }




    /** @brief Return the stencil width associated with the coarsening operator.
     *
     *  The SAMRAI transfer routines guarantee that the source patch will contain
     * sufficient ghostCell data surrounding the interior to satisfy the stencil
     * width requirements for each coarsening operator.
     *
     * In our case, we allow a RF up to 10, so having 5 ghost width is sufficient
     */
    SAMRAI::hier::IntVector getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector{dim, 2};
    }




    /** @brief given a coarseBox, coarse data from the fine patch on the intersection of
     * this box and the box of the destination (the box of the coarse patch).
     *
     * This method will extract fieldData from the two patches, and then
     * get the Field and GridLayout encapsulated into the fieldData.
     * With the help of FieldGeometry, transform the coarseBox to the correct index.
     * After that we can now create FieldCoarsen with the indexAndWeight implementation
     * selected. Finnaly loop over the indexes in the box, and apply the coarsening defined
     * in FieldCoarsen operator
     *
     */
    void coarsen(SAMRAI::hier::Patch& destinationPatch, SAMRAI::hier::Patch const& sourcePatch,
                 int const destinationId, int const sourceId, SAMRAI::hier::Box const& coarseBox,
                 SAMRAI::hier::IntVector const& ratio) const override
    {
        auto& destinationFields  = TensorFieldDataT::getFields(destinationPatch, destinationId);
        auto const& sourceFields = TensorFieldDataT::getFields(sourcePatch, sourceId);
        auto const& sourceLayout = TensorFieldDataT::getLayout(sourcePatch, sourceId);
        auto const& destLayout   = TensorFieldDataT::getLayout(destinationPatch, destinationId);


        // we assume that quantity are the same
        // note that an assertion will be raised in coarseIt operator

        for (std::uint16_t c = 0; c < N; ++c)
        {
            auto const& qty      = destinationFields[c].physicalQuantity();
            using FieldGeometryT = FieldGeometry<GridLayoutT, std::decay_t<decltype(qty)>>;


            // We get different boxes : destination , source, restrictBoxes
            // and transform them in the correct indexing.
            auto const& destPData = destinationPatch.getPatchData(destinationId);
            auto const& srcPData  = sourcePatch.getPatchData(sourceId);
            auto const& destGBox
                = FieldGeometryT::toFieldBox(destPData->getGhostBox(), qty, destLayout);
            auto const& srcGBox
                = FieldGeometryT::toFieldBox(srcPData->getGhostBox(), qty, sourceLayout);
            auto const& coarseLayout   = FieldGeometryT::layoutFromBox(coarseBox, destLayout);
            auto const& coarseFieldBox = FieldGeometryT::toFieldBox(coarseBox, qty, coarseLayout);
            auto const intersectionBox = destGBox * coarseFieldBox;
            // We can now create the coarsening operator
            FieldCoarsenerPolicy coarsener{destLayout.centering(qty), srcGBox, destGBox, ratio};

            coarsen_field(destinationFields[c], sourceFields[c], intersectionBox, coarsener);
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename FieldCoarsenerPolicy,
         typename PhysicalQuantity>
using VecFieldCoarsenOperator
    = TensorFieldCoarsenOperator<1, GridLayoutT, FieldT, FieldCoarsenerPolicy, PhysicalQuantity>;

} // namespace PHARE::amr


#endif
