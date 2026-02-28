#ifndef PHARE_FIELD_REFINE_OPERATOR_HPP
#define PHARE_FIELD_REFINE_OPERATOR_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/box/box_span.hpp"


#include "amr/data/field/field_data.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"


// #include "field_linear_refine.hpp"

#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/hier/RefineOperator.h>


#include <cstddef>


namespace PHARE::amr
{

using core::dirX;
using core::dirY;
using core::dirZ;

template<typename Field_t>
struct CrossLevelIndices
{
    auto constexpr static dim = Field_t::dimension;
    using FieldRow_t          = core::FieldBoxPointRows<Field_t>;
    using FieldRow_ct         = core::FieldBoxPointRows<Field_t const>;


    core::FieldBox<Field_t>& dst;
    core::FieldBox<Field_t const> const& src;

    core::FieldBoxSpan<Field_t, FieldRow_t> dst_lcl_span
        = core::make_field_box_point_span(dst.lcl_box, dst.field);
    core::FieldBoxSpan<Field_t const, FieldRow_ct> src_lcl_span
        = core::make_field_box_point_span(src.lcl_box, src.field);

    core::BoxSpan<int, dim> dst_amr_span = core::make_box_span(dst.amr_box);
    core::BoxSpan<int, dim> src_amr_span = core::make_box_span(src.amr_box);
};




template<typename Field_t>
void refine_field(core::FieldBox<Field_t>& dst, core::FieldBox<Field_t const>& src, auto& refiner)
{
    auto constexpr static IDX = Field_t::dimension - 1;

    CrossLevelIndices<Field_t> indices{dst, src};

    auto d_f_slabs = indices.dst_lcl_span.begin();
    auto d_b_slabs = indices.dst_amr_span.begin();
    auto s_f_slabs = indices.src_lcl_span.begin();
    auto s_b_slabs = indices.src_amr_span.begin();

    auto slab_idx = indices.dst_amr_span.slab_begin();

    for (; d_f_slabs != indices.dst_lcl_span.end(); ++d_f_slabs, ++d_b_slabs, ++slab_idx)
    {
        auto d_f_rows = d_f_slabs.begin();
        auto s_f_rows = s_f_slabs.begin();
        auto d_b_rows = d_b_slabs.begin();
        auto s_b_rows = s_b_slabs.begin();

        auto row_idx = d_b_slabs.row_begin();
        for (; d_f_rows != d_f_slabs.end(); ++d_f_rows, ++d_b_rows, ++row_idx)
        {
            auto const& [d_amr_point, d_size] = *d_b_rows;
            auto const& [s_amr_point, s_size] = *s_b_rows;

            auto&& [d_row, d_lcl_point] = *d_f_rows;
            auto&& [s_row, s_lcl_point] = *s_f_rows;

            std::size_t dst_idx = 0;
            std::size_t src_idx = 0;
            for (; dst_idx < d_row.size(); ++dst_idx)
            {
                assert(s_amr_point == toCoarseIndex(d_amr_point));

                refiner(src.field, dst.field, d_amr_point, s_amr_point, d_lcl_point, s_lcl_point,
                        d_row[dst_idx], s_row[src_idx]);

                if (d_amr_point[IDX] % 2 != 0)
                {
                    ++src_idx;
                    ++s_lcl_point[IDX];
                    ++s_amr_point[IDX];
                }

                ++d_amr_point[IDX];
                ++d_lcl_point[IDX];
            }

            if (row_idx % 2 != 0)
            {
                ++s_f_rows;
                ++s_b_rows;
            }
        }

        if (slab_idx % 2 != 0)
        {
            ++s_f_slabs;
            ++s_b_slabs;
        }
    }
}


template<typename GridLayoutT, typename FieldT, typename FieldRefinerPolicy>
class FieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using GridLayoutImpl                   = typename GridLayoutT::implT;
    using PhysicalQuantity                 = typename FieldT::physical_quantity_type;
    using FieldDataT                       = FieldData<GridLayoutT, FieldT>;

    FieldRefineOperator()
        : SAMRAI::hier::RefineOperator{"FieldRefineOperator"}

    {
    }

    virtual ~FieldRefineOperator() = default;

    /** This implementation have the top priority for refine operation
     *
     */
    NO_DISCARD int getOperatorPriority() const override { return 0; }

    /**
     * @brief This operator needs to have at least 1 ghost cell to work properly
     *
     */
    NO_DISCARD SAMRAI::hier::IntVector
    getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector::getOne(dim);
    }




    /**
     * @brief Given a set of box on a fine patch, compute the interpolation from
     * a coarser patch that is underneath the fine box.
     * Since we get our boxes from a FieldOverlap, we know that they are in correct
     * Field Indexes
     *
     */
    void refine(SAMRAI::hier::Patch& destination, SAMRAI::hier::Patch const& source,
                int const destinationId, int const sourceId,
                SAMRAI::hier::BoxOverlap const& destinationOverlap,
                SAMRAI::hier::IntVector const& ratio) const override
    {
        using FieldGeometry = typename FieldDataT::Geometry;

        auto const& destinationFieldOverlap = dynamic_cast<FieldOverlap const&>(destinationOverlap);

        auto const& overlapBoxes = destinationFieldOverlap.getDestinationBoxContainer();

        auto& destinationField  = FieldDataT::getField(destination, destinationId);
        auto const& destLayout  = FieldDataT::getLayout(destination, destinationId);
        auto const& sourceField = FieldDataT::getField(source, sourceId);
        auto const& srcLayout   = FieldDataT::getLayout(source, sourceId);

        // We assume that quantity are all the same.
        // Note that an assertion will be raised in refineIt operator
        auto const& qty     = destinationField.physicalQuantity();
        auto const destData = destination.getPatchData(destinationId);
        auto const srcData  = source.getPatchData(sourceId);

        auto const destFieldBox
            = FieldGeometry::toFieldBox(destData->getGhostBox(), qty, destLayout);
        auto const sourceFieldBox
            = FieldGeometry::toFieldBox(srcData->getGhostBox(), qty, srcLayout);

        FieldRefinerPolicy refiner{destLayout.centering(qty), destFieldBox, sourceFieldBox, ratio};

        for (auto const& box : overlapBoxes)
        {
            // we compute the intersection with the destination,
            // and then we apply the refine operation on each fine index.
            auto const dst_overlap = phare_box_from<dimension>(destFieldBox * box);
            core::FieldBox dst{destinationField, destLayout, dst_overlap};
            core::FieldBox src{sourceField, srcLayout, coarsen_box(dst_overlap)};
            refine_field(dst, src, refiner);
        }
    }
};


template<std::size_t rank, typename GridLayoutT, typename FieldT, typename FieldRefinerPolicy>
class TensorFieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using GridLayoutImpl                   = GridLayoutT::implT;

    using TensorFieldDataT     = TensorFieldData<rank, GridLayoutT, FieldT, core::HybridQuantity>;
    using TensorFieldOverlap_t = TensorFieldOverlap<rank>;

    static constexpr std::size_t N = TensorFieldDataT::N;

    TensorFieldRefineOperator()
        : SAMRAI::hier::RefineOperator{"TensorFieldRefineOperator"}

    {
    }

    virtual ~TensorFieldRefineOperator() = default;

    /** This implementation have the top priority for refine operation
     *
     */
    NO_DISCARD int getOperatorPriority() const override { return 0; }

    /**
     * @brief This operator needs to have at least 1 ghost cell to work properly
     *
     */
    NO_DISCARD SAMRAI::hier::IntVector
    getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector::getOne(dim);
    }




    /**
     * @brief Given a set of box on a fine patch, compute the interpolation from
     * a coarser patch that is underneath the fine box.
     * Since we get our boxes from a FieldOverlap, we know that they are in correct
     * Field Indexes
     *
     */
    void refine(SAMRAI::hier::Patch& destination, SAMRAI::hier::Patch const& source,
                int const destinationId, int const sourceId,
                SAMRAI::hier::BoxOverlap const& destinationOverlap,
                SAMRAI::hier::IntVector const& ratio) const override
    {
        auto const& destinationTensorFieldOverlap
            = dynamic_cast<TensorFieldOverlap_t const&>(destinationOverlap);
        auto const& srcData      = source.getPatchData(sourceId);
        auto const& destData     = destination.getPatchData(destinationId);
        auto& destinationFields  = TensorFieldDataT::getFields(destination, destinationId);
        auto const& destLayout   = TensorFieldDataT::getLayout(destination, destinationId);
        auto const& sourceFields = TensorFieldDataT::getFields(source, sourceId);
        auto const& srcLayout    = TensorFieldDataT::getLayout(source, sourceId);

        PHARE_LOG_SCOPE(2, "TensorFieldRefineOperator::refine::" + sourceFields[0].name());

        // We assume that quantity are all the same.
        // Note that an assertion will be raised in refineIt operator
        for (std::uint16_t c = 0; c < N; ++c)
        {
            auto const& overlapBoxes
                = destinationTensorFieldOverlap[c]->getDestinationBoxContainer();
            auto const& qty     = destinationFields[c].physicalQuantity();
            using FieldGeometry = FieldGeometry<GridLayoutT, std::decay_t<decltype(qty)>>;

            auto const destFieldBox
                = FieldGeometry::toFieldBox(destData->getGhostBox(), qty, destLayout);
            auto const sourceFieldBox
                = FieldGeometry::toFieldBox(srcData->getGhostBox(), qty, srcLayout);

            FieldRefinerPolicy refiner{destLayout.centering(qty), destFieldBox, sourceFieldBox,
                                       ratio};

            for (auto const& box : overlapBoxes)
            {
                // we compute the intersection with the destination,
                // and then we apply the refine operation on each fine index.
                auto const dst_overlap = phare_box_from<dimension>(destFieldBox * box);
                core::FieldBox dst{destinationFields[c], destLayout, dst_overlap};

                auto const src_overlap
                    = *(coarsen_box(dst_overlap) * srcLayout.AMRGhostBoxFor(sourceFields[c]));
                core::FieldBox src{sourceFields[c], srcLayout, src_overlap};
                refine_field(dst, src, refiner);
            }
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename FieldRefinerPolicy>
using VecFieldRefineOperator
    = TensorFieldRefineOperator</*rank=*/1, GridLayoutT, FieldT, FieldRefinerPolicy>;


} // namespace PHARE::amr



#endif
