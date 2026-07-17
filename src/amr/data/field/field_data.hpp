#ifndef PHARE_SRC_AMR_FIELD_FIELD_DATA_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_DATA_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep


#include "core/logger.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include <core/utilities/types.hpp>
#include "core/utilities/box/box.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/vector.hpp"
#include "core/data/field/field_box.hpp"
#include "core/data/ndarray/ndarray_view.hpp"


#include "amr/samrai.hpp" // restarts
#include <amr/utilities/box/amr_box.hpp>
#include "amr/resources_manager/amr_utils.hpp"


#include "field_geometry.hpp"

#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/tbox/MemoryUtilities.h>

#include <limits>
#include <utility>
#include <stdexcept>



namespace PHARE
{
namespace amr
{

    /**@brief FieldData is the specialization of SAMRAI::hier::PatchData to Field objects
     *
     */
    template<typename GridLayoutT, typename Grid_t,
             typename PhysicalQuantity = decltype(std::declval<Grid_t>().physicalQuantity())>
    class FieldData : public SAMRAI::hier::PatchData
    {
        using Super = SAMRAI::hier::PatchData;

    public:
        using value_type = Grid_t::type;

    private:
        using SetEqualOp          = core::SetEqual<value_type>;
        auto constexpr static NaN = std::numeric_limits<value_type>::quiet_NaN();


        static_assert(core::is_field_v<Grid_t> || core::is_field_tile_set_v<Grid_t>);
        static constexpr auto NO_ROTATE = SAMRAI::hier::Transformation::NO_ROTATE;
        using field_type                = Grid_t::field_type;


    public:
        static constexpr std::size_t dimension    = GridLayoutT::dimension;
        static constexpr std::size_t interp_order = GridLayoutT::interp_order;
        using Geometry                            = FieldGeometry<GridLayoutT, PhysicalQuantity>;
        using gridlayout_type                     = GridLayoutT;
        using data_type                           = field_type;


        /*** \brief Construct a FieldData from information associated to a patch
         *
         * It will create a GridLayout from parameters given by FieldDataFactory
         * From the freshly created GridLayout, it will create a Field with the correct
         * number of cells in each needed directions
         */
        FieldData(SAMRAI::hier::Box const& domain, SAMRAI::hier::IntVector const& ghost,
                  std::string name, GridLayoutT const& layout, PhysicalQuantity qty)
            : SAMRAI::hier::PatchData(domain, ghost)
            , gridLayout{layout}
            , field(name, layout, qty, NaN)
            , quantity_{qty}
        {
        } //


        FieldData()                            = delete;
        FieldData(FieldData const&)            = delete;
        FieldData(FieldData&&)                 = default;
        FieldData& operator=(FieldData const&) = delete;



        void getFromRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) override
        {
            Super::getFromRestart(restart_db);

            if constexpr (core::is_field_tile_set_v<Grid_t>)
                throw std::runtime_error("finish");
            else
                getVectorFromRestart(*restart_db, "field_" + field.name(), field.vector());
        }

        void putToRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) const override
        {
            Super::putToRestart(restart_db);

            if constexpr (core::is_field_tile_set_v<Grid_t>)
                throw std::runtime_error("finish");
            else
                putVectorToRestart(*restart_db, "field_" + field.name(), field.vector());
        }



        /*** \brief Copy information from another FieldData where data overlap
         *
         *    The data will be copied from the interior and ghost of the source to the interior
         * and ghost of the destination, where there is an overlap in the underlying index space
         */
        void copy(SAMRAI::hier::PatchData const& source) final
        {
            PHARE_LOG_SCOPE(3, "FieldData::copy");

            // After checking that source and *this have the same number of dimension
            // We will try to cast source as a FieldData, if it succeed we can continue
            // and perform the copy. Otherwise we call copy2 that will simply throw a runtime
            // error

            TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source);

            // throws on failure
            auto& fieldSource = dynamic_cast<FieldData const&>(source);

            TBOX_ASSERT(quantity_ == fieldSource.quantity_);
            // First step is to translate the AMR box into proper index space of the given
            // quantity_ using the source gridlayout to accomplish that we get the interior box,
            // from the FieldData.

            SAMRAI::hier::Box const sourceBox = Geometry::toFieldBox(
                fieldSource.getGhostBox(), quantity_, fieldSource.gridLayout);

            SAMRAI::hier::Box const destinationBox
                = Geometry::toFieldBox(this->getGhostBox(), quantity_, this->gridLayout);

            // Given the two boxes in correct space we just have to intersect them
            SAMRAI::hier::Box const intersectionBox = sourceBox * destinationBox;

            if (!intersectionBox.empty())
            {
                // We can copy field from the source to the destination on the correct region
                copy_(intersectionBox, sourceBox, destinationBox, fieldSource, field);
            }
        }




        /*** \brief This form should not be called since we cannot derive from FieldData
         * since FieldData is a final implementation of PatchData
         */
        void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination) const final
        {
            throw std::runtime_error("Error cannot cast the PatchData to FieldData");
        }




        /*** \brief Copy data from the source into the destination using the designated overlap
         * descriptor.
         *
         *   The overlap will contain AMR index space boxes on destination to be filled and also
         * give the necessary transformation to apply to the source, to perform the copy (ie :
         * translation for periodics condition)
         */
        void copy(SAMRAI::hier::PatchData const& source,
                  SAMRAI::hier::BoxOverlap const& overlap) final
        {
            PHARE_LOG_SCOPE(3, "FieldData::copy");

            // casts throw on failure
            auto& fieldSource  = dynamic_cast<FieldData const&>(source);
            auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

            copy_(fieldSource, fieldOverlap);
        }




        /*** \brief This form should not be called since we cannot derive from FieldData
         */
        void copy2([[maybe_unused]] SAMRAI::hier::PatchData& destination,
                   [[maybe_unused]] SAMRAI::hier::BoxOverlap const& overlap) const final
        {
            throw std::runtime_error("Error cannot cast the PatchData to FieldData");
        }




        /*** \brief Determines whether the patch data subclass can estimate the necessary stream
         * size using only index space information.
         *
         * The return value is true since that for a corresponding domain, there is a fixed
         * number of elements in the field depending on the PhysicalQuantity and the Layout used
         */
        bool canEstimateStreamSizeFromBox() const final { return true; }



        /*** \brief Compute the maximum amount of memory needed to hold FieldData information on
         * the specified overlap
         */
        std::size_t getDataStreamSize(SAMRAI::hier::BoxOverlap const& overlap) const final
        {
            return getDataStreamSize_(overlap);
        }




        /*** \brief Serialize the data contained in the field data on the region covered by the
         * overlap, and put it on the stream.
         */
        void packStream(SAMRAI::tbox::MessageStream& stream,
                        SAMRAI::hier::BoxOverlap const& overlap) const final
        {
            packStream(stream, overlap, field);
        }

        template<typename Operator = SetEqualOp>
        void packStream(SAMRAI::tbox::MessageStream& stream,
                        SAMRAI::hier::BoxOverlap const& overlap, auto const& src) const
        {
            PHARE_LOG_SCOPE(3, "packStream");

            auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

            SAMRAI::hier::Transformation const& transformation = fieldOverlap.getTransformation();

            if (transformation.getRotation() != NO_ROTATE)
                throw std::runtime_error("Rotations are not supported in PHARE");

            auto const buffer_size = getDataStreamSize_(overlap) / sizeof(value_type);
            auto& buffer           = tmp.reserve_and_clear(buffer_size)();

            for (auto const& box : fieldOverlap.getDestinationBoxContainer())
            {
                SAMRAI::hier::Box packBox{box};

                // Since the transformation, allow to transform the source box,
                // into the destination box space, and that the box in the boxContainer
                // are in destination space, we have to use the inverseTransform
                // to get into source space
                transformation.inverseTransform(packBox);

                core::FieldBox<Grid_t const>{src, gridLayout, phare_box_from<dimension>(packBox)}
                    .template append_to<Operator>(buffer);
            }

            // Once we have fill the buffer, we send it on the stream
            stream.pack(buffer.data(), buffer.size());
        }



        /*** \brief Unserialize data contained on the stream, that comes from a region covered
         * by the overlap, and fill the data where is needed.
         */
        void unpackStream(SAMRAI::tbox::MessageStream& stream,
                          SAMRAI::hier::BoxOverlap const& overlap) final
        {
            unpackStream(stream, overlap, field);
        }

        template<typename Operator = SetEqualOp>
        void unpackStream(SAMRAI::tbox::MessageStream& stream,
                          SAMRAI::hier::BoxOverlap const& overlap, Grid_t& dst)
        {
            PHARE_LOG_SCOPE(3, "unpackStream");

            auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

            if (fieldOverlap.getTransformation().getRotation() != NO_ROTATE)
                throw std::runtime_error("Rotations are not supported in PHARE");

            // For unpacking we need to know how much element we will need to extract
            auto& buffer = tmp.get_no_copy(getDataStreamSize(overlap) / sizeof(value_type))();

            // We flush a portion of the stream on the buffer.
            stream.unpack(buffer.data(), buffer.size());

            // Here the seek counter will be used to index buffer
            std::size_t seek = 0;

            // For unpackStream, there is no transformation needed, since all the box
            // are on the destination space
            for (auto const& sambox : fieldOverlap.getDestinationBoxContainer())
            {
                auto const box = phare_box_from<dimension>(sambox);
                auto const view
                    = core::make_array_view(buffer.data() + seek, *box.shape().as_unsigned());
                using NdArray = decltype(view);

                core::FieldBox<Grid_t>{dst, gridLayout, box}.template op<Operator>(
                    core::FieldBox<NdArray>{
                        view, box,
                        core::box_from_zero_to_upper_minus_one(*box.shape().as_unsigned())});
                seek += sambox.size();
            }
        }



        auto* getPointer() { return &field; }


        static GridLayoutT const& getLayout(SAMRAI::hier::Patch const& patch, int id)
        {
            auto const& patchData
                = std::dynamic_pointer_cast<FieldData<GridLayoutT, Grid_t>>(patch.getPatchData(id));
            if (!patchData)
            {
                throw std::runtime_error("cannot cast to FieldData");
            }
            return patchData->gridLayout;
        }


        static Grid_t& getField(SAMRAI::hier::Patch const& patch, int id)
        {
            auto const& patchData
                = std::dynamic_pointer_cast<FieldData<GridLayoutT, Grid_t>>(patch.getPatchData(id));
            if (!patchData)
            {
                throw std::runtime_error("cannot cast to FieldData");
            }
            return patchData->field;
        }


        template<typename Operation>
        void operate(auto const& src, SAMRAI::hier::BoxOverlap const& overlap);


        template<typename Operation>
        void packStreamAnd(SAMRAI::tbox::MessageStream& stream,
                           SAMRAI::hier::BoxOverlap const& overlap);

        template<typename Operation>
        void unpackStreamAnd(SAMRAI::tbox::MessageStream& stream,
                             SAMRAI::hier::BoxOverlap const& overlap);


        GridLayoutT gridLayout;
        Grid_t field;

    private:
        PhysicalQuantity quantity_; ///! PhysicalQuantity used for this field data
        static inline core::MinimizingVector<value_type> tmp; // LESS ALLOCATIONS




        /*** \brief copy data from the intersection box
         *
         */

        template<typename Operator = SetEqualOp>
        void copy_(SAMRAI::hier::Box const& intersectBox, SAMRAI::hier::Box const& sourceBox,
                   SAMRAI::hier::Box const& destinationBox, FieldData const& source,
                   Grid_t& fieldDestination)
        {
            // First we represent the intersection that is defined in AMR space to the local
            // space of the source Then we represent the intersection into the local space of
            // the destination We can finally perform the copy of the element in the correct
            // range

            core::FieldBox<Grid_t>{
                fieldDestination, gridLayout,
                as_unsigned_phare_box<dimension>(AMRToLocal(intersectBox, destinationBox))}
                .template op<Operator>(core::FieldBox<Grid_t const>{
                    source.field, source.gridLayout,
                    as_unsigned_phare_box<dimension>(AMRToLocal(intersectBox, sourceBox))});
        }


        void copy_(FieldData const& source, FieldOverlap const& overlap)
        {
            copy_(source, overlap, field);
        }



        template<typename Operator = SetEqualOp, typename... T0>
        void copy_(FieldData<T0...> const& source, FieldOverlap const& overlap, Grid_t& dst)
        {
            // Here the first step is to get the transformation from the overlap
            // we transform the box from the source, and from the destination
            // from AMR index to FieldData indexes (ie whether or not the quantity is primal
            // or not), and we also consider the ghost. After that we compute the
            // intersection with the source box, the destinationBox, and the box from the
            // destinationBoxContainer.

            SAMRAI::hier::Transformation const& transformation = overlap.getTransformation();

            if (transformation.getRotation() != NO_ROTATE)
                throw std::runtime_error("copy with rotate not implemented");

            if (transformation.getBeginBlock() == transformation.getEndBlock())
            {
                auto const& offset = as_point<dimension>(transformation);

                for (auto const& box : overlap.getDestinationBoxContainer())
                {
                    SAMRAI::hier::Box const sourceBox
                        = Geometry::toFieldBox(source.getGhostBox(), quantity_, source.gridLayout);

                    SAMRAI::hier::Box const destinationBox
                        = Geometry::toFieldBox(this->getGhostBox(), quantity_, this->gridLayout);

                    SAMRAI::hier::Box transformedSource{sourceBox};
                    transformation.transform(transformedSource);

                    SAMRAI::hier::Box const intersectionBox{box * transformedSource
                                                            * destinationBox};

                    if (!intersectionBox.empty())
                        core::FieldBox{dst, gridLayout,
                                       as_unsigned_phare_box<dimension>(
                                           AMRToLocal(intersectionBox, destinationBox))}
                            .template op<Operator>(
                                core::FieldBox{source.field, source.gridLayout,
                                               as_unsigned_phare_box<dimension>(
                                                   AMRToLocal(intersectionBox, transformedSource))}
                                    .offset(offset));
                }
            }
        }


        /*** \brief Compute the maximum amount of memory needed to hold FieldData information on
         *          the specified overlap
         */
        std::size_t getDataStreamSize_(SAMRAI::hier::BoxOverlap const& overlap) const
        {
            // The idea here is to tell SAMRAI the maximum memory will be used by our type
            // on a given region.

            // throws on failure
            auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

            if (fieldOverlap.isOverlapEmpty())
                return 0;

            // TODO: see FieldDataFactory todo of the same function

            SAMRAI::hier::BoxContainer const& boxContainer
                = fieldOverlap.getDestinationBoxContainer();

            return boxContainer.getTotalSizeOfBoxes() * sizeof(value_type);
        }
    };



} // namespace amr
} // namespace PHARE



namespace PHARE::amr
{

template<typename GridLayoutT, typename Grid_t, typename PhysicalQuantity>
template<typename Operation>
void FieldData<GridLayoutT, Grid_t, PhysicalQuantity>::packStreamAnd(
    SAMRAI::tbox::MessageStream& stream, SAMRAI::hier::BoxOverlap const& overlap)
{
    packStream<Operation>(stream, overlap, field);
}

template<typename GridLayoutT, typename Grid_t, typename PhysicalQuantity>
template<typename Operation>
void FieldData<GridLayoutT, Grid_t, PhysicalQuantity>::unpackStreamAnd(
    SAMRAI::tbox::MessageStream& stream, SAMRAI::hier::BoxOverlap const& overlap)
{
    unpackStream<Operation>(stream, overlap, field);
}

template<typename GridLayoutT, typename Grid_t, typename PhysicalQuantity>
template<typename Operation>
void FieldData<GridLayoutT, Grid_t, PhysicalQuantity>::operate(
    auto const& src, SAMRAI::hier::BoxOverlap const& overlap)
{
    TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

    auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

    copy_<Operation>(src, fieldOverlap, field);
}


} // namespace PHARE::amr


#endif
