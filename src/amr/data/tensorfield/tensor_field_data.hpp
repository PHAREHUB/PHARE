#ifndef PHARE_SRC_AMR_TENSORFIELD_TENSORFIELD_DATA_HPP
#define PHARE_SRC_AMR_TENSORFIELD_TENSORFIELD_DATA_HPP

#include "amr/resources_manager/amr_utils.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/logger.hpp"
#include "core/data/field/field_box.hpp"
#include "core/data/tensorfield/tensorfield.hpp"


#include "amr/data/field/field_overlap.hpp"
#include "amr/data/tensorfield/tensor_field_geometry.hpp"

#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/tbox/MemoryUtilities.h>
#include <optional>


namespace PHARE::amr
{
// We use another class here so that we can specialize specifics function: copy , pack , unpack
// on the dimension and we don't want to loose non specialized function related to SAMRAI
// interface
template<typename GridLayoutT, std::size_t dim, typename Grid_t, typename PhysicalQuantity>
class TensorFieldDataInternals
{
};

/**
 * @brief TensorFieldData is the specialization of SAMRAI::hier::PatchData to Field objects
 */
template<std::size_t rank, typename GridLayoutT, typename Grid_t, typename PhysicalQuantity>
class TensorFieldData : public SAMRAI::hier::PatchData
{
    using This  = TensorFieldData<rank, GridLayoutT, Grid_t, PhysicalQuantity>;
    using Super = SAMRAI::hier::PatchData;

    static constexpr auto NO_ROTATE = SAMRAI::hier::Transformation::NO_ROTATE;

    using tensor_t = typename PhysicalQuantity::template TensorType<rank>;

    template<typename ComponentNames, typename GridLayout>
    auto static make_grids(ComponentNames const& compNames, GridLayout const& layout, tensor_t qty)
    {
        auto qts = PhysicalQuantity::componentsQuantities(qty);
        return core::for_N<N, core::for_N_R_mode::make_array>(
            [&](auto i) { return Grid_t{compNames[i], qts[i], layout.allocSize(qts[i])}; });
    }

    using value_type = Grid_t::value_type;
    using SetEqualOp = core::Equals<value_type>;

    template<typename Box_t>
    std::optional<Box_t> dual_dir_minus_1(Box_t box, auto const qty) const
    {
        auto const centerings = GridLayoutT::centering(qty);
        for (std::uint16_t i = 0; i < dimension; ++i)
        {
            auto const dual_dir = centerings[i] == core::QtyCentering::dual;
            if (dual_dir)
            {
                if (box.lower[i] == box.upper[i]) // bad!
                    return std::nullopt;

                box.upper[i] -= 1;
            }
        }
        return box;
    }

public:
    static constexpr std::size_t dimension    = GridLayoutT::dimension;
    static constexpr std::size_t interp_order = GridLayoutT::interp_order;
    static constexpr auto N                   = core::detail::tensor_field_dim_from_rank<rank>();

    using Geometry        = TensorFieldGeometry<rank, GridLayoutT, PhysicalQuantity>;
    using gridlayout_type = GridLayoutT;

    /*** \brief Construct a TensorFieldData from information associated to a patch
     *
     * It will create a GridLayout from parameters given by TensorFieldDataFactory
     * From the freshly created GridLayout, it will create a Field with the correct
     * number of cells in each needed directions
     */
    TensorFieldData(SAMRAI::hier::Box const& domain, SAMRAI::hier::IntVector const& ghost,
                    std::string name, GridLayoutT const& layout, tensor_t qty)
        : SAMRAI::hier::PatchData(domain, ghost)
        , gridLayout{layout}
        , grids(make_grids(core::detail::tensor_field_names<rank>(name), layout, qty))
        , quantity_{qty}
    {
    }


    TensorFieldData()                                  = delete;
    TensorFieldData(TensorFieldData const&)            = delete;
    TensorFieldData(TensorFieldData&&)                 = default;
    TensorFieldData& operator=(TensorFieldData const&) = delete;



    void getFromRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) override
    {
        Super::getFromRestart(restart_db);

        for (std::uint16_t c = 0; c < N; ++c)
        {
            assert(grids[c].vector().size() > 0);
            restart_db->getDoubleArray("field_" + grids[c].name(), grids[c].vector().data(),
                                       grids[c].vector().size()); // do not reallocate!
        }
    }

    void putToRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) const override
    {
        Super::putToRestart(restart_db);

        for (std::uint16_t c = 0; c < N; ++c)
            restart_db->putVector("field_" + grids[c].name(), grids[c].vector());
    };




    /*** \brief Copy information from another TensorFieldData where data overlap
     *
     *    The data will be copied from the interior and ghost of the source to the interior and
     *    ghost of the destination, where there is an overlap in the underlying index space
     */
    void copy(const SAMRAI::hier::PatchData& source) final
    {
        PHARE_LOG_SCOPE(3, "TensorFieldData::copy");

        // After checking that source and *this have the same number of dimension
        // We will try to cast source as a TensorFieldData, if it succeed we can continue
        // and perform the copy. Otherwise we call copy2 that will simply throw a runtime
        // error

        TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source);

        // throws on failure
        auto& fieldSource = dynamic_cast<TensorFieldData const&>(source);

        TBOX_ASSERT(quantity_ == fieldSource.quantity_);
        // First step is to translate the AMR box into proper index space of the given
        // quantity_ using the source gridlayout to accomplish that we get the interior box,
        // from the TensorFieldData.

        SAMRAI::hier::Box sourceBox
            = Geometry::toFieldBox(fieldSource.getGhostBox(), quantity_, fieldSource.gridLayout);


        SAMRAI::hier::Box destinationBox
            = Geometry::toFieldBox(this->getGhostBox(), quantity_, this->gridLayout);

        // Given the two boxes in correct space we just have to intersect them
        SAMRAI::hier::Box intersectionBox = sourceBox * destinationBox;

        if (!intersectionBox.empty())
            copy_(intersectionBox, sourceBox, destinationBox, fieldSource, *this);
    }




    /*** \brief This form should not be called since we cannot derive from TensorFieldData
     * since TensorFieldData is a final implementation of PatchData
     */
    void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination) const final
    {
        throw std::runtime_error("Error cannot cast the PatchData to TensorFieldData");
    }




    /*** \brief Copy data from the source into the destination using the designated overlap
     * descriptor.
     *
     *   The overlap will contain AMR index space boxes on destination to be filled and also
     * give the necessary transformation to apply to the source, to perform the copy (ie :
     * translation for periodics condition)
     */
    void copy(const SAMRAI::hier::PatchData& source, const SAMRAI::hier::BoxOverlap& overlap) final
    {
        PHARE_LOG_SCOPE(3, "TensorFieldData::copy");

        // casts throw on failure
        auto& fieldSource  = dynamic_cast<TensorFieldData const&>(source);
        auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

        copy_(fieldSource, fieldOverlap);
    }




    /*** \brief This form should not be called since we cannot derive from TensorFieldData
     */
    void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination,
               [[maybe_unused]] const SAMRAI::hier::BoxOverlap& overlap) const final
    {
        throw std::runtime_error("Error cannot cast the PatchData to TensorFieldData");
    }




    /*** \brief Determines whether the patch data subclass can estimate the necessary stream
     * size using only index space information.
     *
     * The return value is true since that for a corresponding domain, there is a fixed
     * number of elements in the field depending on the PhysicalQuantity and the Layout used
     */
    bool canEstimateStreamSizeFromBox() const final { return true; }



    /*** \brief Compute the maximum amount of memory needed to hold TensorFieldData information
     * on the specified overlap
     */
    std::size_t getDataStreamSize(const SAMRAI::hier::BoxOverlap& overlap) const final
    {
        return getDataStreamSize_(overlap);
    }




    /*** \brief Serialize the data contained in the field data on the region covered by the
     * overlap, and put it on the stream.
     */
    void packStream(SAMRAI::tbox::MessageStream& stream,
                    const SAMRAI::hier::BoxOverlap& overlap) const final
    {
        PHARE_LOG_SCOPE(3, "packStream");

        std::size_t const expectedSize = getDataStreamSize_(overlap) / sizeof(double);
        std::vector<typename Grid_t::type> buffer;
        buffer.reserve(expectedSize * N); // :(

        auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

        SAMRAI::hier::Transformation const& transformation = fieldOverlap.getTransformation();
        if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
        {
            for (auto const& box : fieldOverlap.getDestinationBoxContainer())
            {
                for (std::size_t c = 0; c < N; ++c)
                {
                    auto const& source = grids[c];
                    SAMRAI::hier::Box packBox{box};

                    // Since the transformation, allow to transform the source box,
                    // into the destination box space, and that the box in the boxContainer
                    // are in destination space, we have to use the inverseTransform
                    // to get into source space
                    transformation.inverseTransform(packBox);

                    if (auto const finalBox = dual_dir_minus_1(phare_box_from<dimension>(packBox),
                                                               source.physicalQuantity()))
                    {
                        core::FieldBox<Grid_t const> src{source, gridLayout, *finalBox};
                        src.append_to(buffer);
                    }
                }
            }
        }
        // throw, we don't do rotations in phare....

        // Once we have fill the buffer, we send it on the stream
        stream.pack(buffer.data(), buffer.size());
    }




    /*** \brief Unserialize data contained on the stream, that comes from a region covered by
     * the overlap, and fill the data where is needed.
     */
    void unpackStream(SAMRAI::tbox::MessageStream& stream,
                      const SAMRAI::hier::BoxOverlap& overlap) final
    {
        unpackStream(stream, overlap, grids);
    }

    template<typename Operator = SetEqualOp>
    void unpackStream(SAMRAI::tbox::MessageStream& stream, const SAMRAI::hier::BoxOverlap& overlap,
                      auto& dst_grids)
    {
        PHARE_LOG_SCOPE(3, "unpackStream");

        auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

        if (fieldOverlap.getTransformation().getRotation() != NO_ROTATE)
            throw std::runtime_error("Rotations are not supported in PHARE");

        // For unpacking we need to know how much element we will need to extract
        std::vector<double> buffer(getDataStreamSize(overlap) / sizeof(value_type), 0.);

        // We flush a portion of the stream on the buffer.
        stream.unpack(buffer.data(), buffer.size());

        // Here the seek counter will be used to index buffer
        std::size_t seek = 0;

        // For unpackStream, there is no transformation needed, since all the box
        // are on the destination space

        for (std::size_t c = 0; c < N; ++c)
            for (auto const& sambox : fieldOverlap.getDestinationBoxContainer())
            {
                auto& dst_grid = dst_grids[c];
                if (auto const box = dual_dir_minus_1(phare_box_from<dimension>(sambox),
                                                      dst_grid.physicalQuantity()))
                {
                    core::FieldBox<Grid_t> dst{dst_grid, gridLayout, *box};
                    dst.template set_from<Operator>(buffer, seek);
                    seek += box->size();
                }
            }
    }



    auto* getPointer() { return &grids; }


    static GridLayoutT const& getLayout(SAMRAI::hier::Patch const& patch, int id)
    {
        auto const& patchData = std::dynamic_pointer_cast<This>(patch.getPatchData(id));
        if (!patchData)
            throw std::runtime_error("cannot cast to TensorFieldData");
        return patchData->gridLayout;
    }


    static auto& getFields(SAMRAI::hier::Patch const& patch, int const id)
    {
        auto const& patchData = std::dynamic_pointer_cast<This>(patch.getPatchData(id));
        if (!patchData)
            throw std::runtime_error("cannot cast to TensorFieldData");
        return patchData->grids;
    }

    void sum(SAMRAI::hier::PatchData const& src, SAMRAI::hier::BoxOverlap const& overlap);
    void unpackStreamAndSum(SAMRAI::tbox::MessageStream& stream,
                            SAMRAI::hier::BoxOverlap const& overlap);



    GridLayoutT gridLayout;
    std::array<Grid_t, N> grids;

private:
    tensor_t quantity_; ///! PhysicalQuantity used for this field data




    /*** \brief copy data from the intersection box
     *
     */
    template<typename Operator = SetEqualOp>
    void copy_(SAMRAI::hier::Box const& intersectBox, SAMRAI::hier::Box const& sourceBox,
               SAMRAI::hier::Box const& destinationBox, TensorFieldData const& source,
               TensorFieldData& destination)
    {
        // First we represent the intersection that is defined in AMR space to the local
        // space of the source Then we represent the intersection into the local space of
        // the destination We can finally perform the copy of the element in the correct
        // range

        for (std::size_t c = 0; c < N; ++c)
        {
            auto& dst_grid     = destination.grids[c];
            auto& src_grid     = source.grids[c];
            auto const dst_box = dual_dir_minus_1(
                as_unsigned_phare_box<dimension>(AMRToLocal(intersectBox, destinationBox)),
                dst_grid.physicalQuantity());
            auto const src_box = dual_dir_minus_1(
                as_unsigned_phare_box<dimension>(AMRToLocal(intersectBox, sourceBox)),
                src_grid.physicalQuantity());
            if (dst_box and src_box)
            {
                core::FieldBox<Grid_t> dst{dst_grid, gridLayout, *dst_box};
                core::FieldBox<Grid_t const> const src{src_grid, source.gridLayout, *src_box};
                operate_on_fields<Operator>(dst, src);
            }
        }
    }


    void copy_(TensorFieldData const& source, FieldOverlap const& overlap)
    {
        copy_(source, overlap, *this);
    }

    template<typename Operator = SetEqualOp>
    void copy_(TensorFieldData const& source, FieldOverlap const& overlap, TensorFieldData& dst)
    {
        // Here the first step is to get the transformation from the overlap
        // we transform the box from the source, and from the destination
        // from AMR index to TensorFieldData indexes (ie whether or not the quantity is primal
        // or not), and we also consider the ghost. After that we compute the
        // intersection with the source box, the destinationBox, and the box from the
        // destinationBoxContainer.


        SAMRAI::hier::Transformation const& transformation = overlap.getTransformation();

        if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
        {
            SAMRAI::hier::BoxContainer const& boxList = overlap.getDestinationBoxContainer();

            SAMRAI::hier::IntVector const zeroOffset{
                SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension{dimension})};

            if (transformation.getBeginBlock() == transformation.getEndBlock())
            {
                for (auto const& box : boxList)
                {
                    SAMRAI::hier::Box sourceBox
                        = Geometry::toFieldBox(source.getGhostBox(), quantity_, source.gridLayout);


                    SAMRAI::hier::Box destinationBox
                        = Geometry::toFieldBox(this->getGhostBox(), quantity_, this->gridLayout);


                    SAMRAI::hier::Box transformedSource{sourceBox};
                    transformation.transform(transformedSource);


                    SAMRAI::hier::Box intersectionBox{box * transformedSource * destinationBox};


                    if (!intersectionBox.empty())
                        copy_<Operator>(intersectionBox, transformedSource, destinationBox, source,
                                        dst);
                }
            }
        }
        else
        {
            throw std::runtime_error("copy with rotate not implemented");
        }
    }



    std::size_t getDataStreamSize_(SAMRAI::hier::BoxOverlap const& overlap) const
    {
        // The idea here is to tell SAMRAI the maximum memory will be used by our type
        // on a given region.


        // throws on failure
        auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

        if (fieldOverlap.isOverlapEmpty())
        {
            return 0;
        }

        SAMRAI::hier::BoxContainer const& boxContainer = fieldOverlap.getDestinationBoxContainer();

        return boxContainer.getTotalSizeOfBoxes() * sizeof(typename Grid_t::type) * N;
    }


}; // namespace PHARE




template<std::size_t rank, typename GridLayoutT, typename Grid_t, typename PhysicalQuantity>
void TensorFieldData<rank, GridLayoutT, Grid_t, PhysicalQuantity>::unpackStreamAndSum(
    SAMRAI::tbox::MessageStream& stream, SAMRAI::hier::BoxOverlap const& overlap)
{
    using PlusEqualOp = core::PlusEquals<value_type>;

    unpackStream<PlusEqualOp>(stream, overlap, grids);
}



template<std::size_t rank, typename GridLayoutT, typename Grid_t, typename PhysicalQuantity>
void TensorFieldData<rank, GridLayoutT, Grid_t, PhysicalQuantity>::sum(
    SAMRAI::hier::PatchData const& src, SAMRAI::hier::BoxOverlap const& overlap)
{
    using PlusEqualOp = core::PlusEquals<value_type>;

    TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

    auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);
    auto& fieldSource  = dynamic_cast<TensorFieldData const&>(src);

    copy_<PlusEqualOp>(fieldSource, fieldOverlap, *this);
}


} // namespace PHARE::amr


#endif
