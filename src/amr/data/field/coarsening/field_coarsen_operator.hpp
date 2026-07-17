#ifndef PHARE_FIELD_DATA_COARSEN_HPP
#define PHARE_FIELD_DATA_COARSEN_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "amr/data/field/field_data.hpp"
#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/CoarsenOperator.h>


namespace PHARE::amr
{


template<typename Coarsener, typename FieldT>
void coarsen_field(FieldT& dst, auto const& dstGhostBox, auto const& dstLayout, FieldT const& src,
                   auto const& srcGhostBox, auto& box, auto ratio)
    requires(not core::is_field_tile_set_v<FieldT>)
{
    auto const& qty = dst.physicalQuantity();

    Coarsener coarsener{dstLayout.centering(qty), srcGhostBox, dstGhostBox, ratio};
    for (auto const bix : box)
        coarsener(src, dst, bix);
}



template<typename Coarsener, typename FieldT>
void coarsen_field(FieldT& dst, auto const& dstGhostBox, auto const& dstLayout, FieldT const& src,
                   auto const& srcGhostBox, auto& box, auto ratio)
    requires(core::is_field_tile_set_v<FieldT>)
{
    auto const& qty = dst.physicalQuantity();
    for (auto& dst_tile : dst())
    {
        auto const dst_box = dst_tile.ghost_box();
        if (auto const dst_overlap = dst_box * box)
        {
            for (auto const& src_tile : src())
            {
                auto const src_box = src_tile.field_box();
                if (auto const src_overlap = src_box * refine_box(box))
                {
                    if (auto const& overlap = *dst_overlap * coarsen_box(*src_overlap))
                    {
                        Coarsener coarsener{dstLayout.centering(qty),
                                            samrai_box_from(src_tile.ghost_box()),
                                            samrai_box_from(dst_tile.ghost_box()), ratio};

                        for (auto const bix : *overlap)
                            coarsener(src_tile(), dst_tile(), bix);
                    }
                }
            }
        }
    }
}



} // namespace PHARE::amr


namespace PHARE
{
namespace amr
{

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
        FieldCoarsenOperator& operator=(FieldCoarsenOperator&&)      = delete;


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
         * selected. Finally loop over the indexes in the box, and apply the coarsening defined
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

            auto const box = phare_box_from<dimension>(intersectionBox);

            coarsen_field<FieldCoarsenerPolicy>( //
                destinationField, destGBox, destLayout, sourceField, srcGBox, box, ratio);
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

            auto const box = phare_box_from<dimension>(intersectionBox);
            coarsen_field<FieldCoarsenerPolicy>( //
                destinationFields[c], destGBox, destLayout, sourceFields[c], srcGBox, box, ratio);
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename FieldCoarsenerPolicy,
         typename PhysicalQuantity>
using VecFieldCoarsenOperator = TensorFieldCoarsenOperator</*rank=*/1, GridLayoutT, FieldT,
                                                           FieldCoarsenerPolicy, PhysicalQuantity>;

} // namespace PHARE::amr


#endif
