#ifndef PHARE_SRC_AMR_FIELD_FIELD_DATA_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_DATA_HPP



#include "core/def/phare_mpi.hpp"

#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/tbox/MemoryUtilities.h>
#include <utility>

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "field_geometry.hpp"

#include "core/logger.hpp"

#include <iostream>

namespace PHARE
{
namespace amr
{
    // We use another class here so that we can specialize specifics function: copy , pack , unpack
    // on the dimension and we don't want to loose non specialized function related to SAMRAI
    // interface
    template<typename GridLayoutT, std::size_t dim, typename FieldImpl,
             typename PhysicalQuantity = decltype(std::declval<FieldImpl>().physicalQuantity())>
    class FieldDataInternals
    {
    };


    /**@brief FieldData is the specialization of SAMRAI::hier::PatchData to Field objects
     *
     */
    template<typename GridLayoutT, typename FieldImpl,
             typename PhysicalQuantity = decltype(std::declval<FieldImpl>().physicalQuantity())>
    class FieldData : public SAMRAI::hier::PatchData
    {
        using Super = SAMRAI::hier::PatchData;

    public:
        static constexpr std::size_t dimension    = GridLayoutT::dimension;
        static constexpr std::size_t interp_order = GridLayoutT::interp_order;
        using Geometry                            = FieldGeometry<GridLayoutT, PhysicalQuantity>;

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
            , field(name, qty, gridLayout.allocSize(qty))
            , quantity_{qty}
        {
        } //

        [[deprecated]] FieldData(SAMRAI::hier::Box const& domain,
                                 SAMRAI::hier::IntVector const& ghost, std::string name,
                                 std::array<double, dimension> const& dl,
                                 std::array<std::uint32_t, dimension> const& nbrCells,
                                 core::Point<double, dimension> const& origin, PhysicalQuantity qty)

            : SAMRAI::hier::PatchData(domain, ghost)
            , gridLayout{dl, nbrCells, origin}
            , field(name, qty, gridLayout.allocSize(qty))
            , quantity_{qty}
        {
        }

        FieldData()                            = delete;
        FieldData(FieldData const&)            = delete;
        FieldData(FieldData&&)                 = default;
        FieldData& operator=(FieldData const&) = delete;



        void getFromRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) override
        {
            Super::getFromRestart(restart_db);

            restart_db->getVector("field_" + field.name(), field.vector());
        }

        void putToRestart(std::shared_ptr<SAMRAI::tbox::Database> const& restart_db) const override
        {
            Super::putToRestart(restart_db);

            restart_db->putVector("field_" + field.name(), field.vector());
        };




        /*** \brief Copy information from another FieldData where data overlap
         *
         *    The data will be copied from the interior and ghost of the source to the interior and
         *    ghost of the destination, where there is an overlap in the underlying index space
         */
        void copy(const SAMRAI::hier::PatchData& source) final
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
            // from the FieldData. and we call toFieldBox (with the default parameter withGhost
            // = true). note that we could have stored the ghost box of the field data at
            // creation

            SAMRAI::hier::Box sourceBox
                = Geometry::toFieldBox(fieldSource.getBox(), quantity_, fieldSource.gridLayout);


            SAMRAI::hier::Box destinationBox
                = Geometry::toFieldBox(this->getBox(), quantity_, this->gridLayout);

            // Given the two boxes in correct space we just have to intersect them
            SAMRAI::hier::Box intersectionBox = sourceBox * destinationBox;

            if (!intersectionBox.empty())
            {
                auto const& sourceField = fieldSource.field;
                auto& destinationField  = field;

                // We can copy field from the source to the destination on the correct region
                copy_(intersectionBox, sourceBox, destinationBox, fieldSource, sourceField,
                      destinationField);
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
        void copy(const SAMRAI::hier::PatchData& source,
                  const SAMRAI::hier::BoxOverlap& overlap) final
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
                   [[maybe_unused]] const SAMRAI::hier::BoxOverlap& overlap) const final
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
        std::size_t getDataStreamSize(const SAMRAI::hier::BoxOverlap& overlap) const final
        {
            return getDataStreamSize_<false>(overlap);
        }




        /*** \brief Serialize the data contained in the field data on the region covered by the
         * overlap, and put it on the stream.
         */
        void packStream(SAMRAI::tbox::MessageStream& stream,
                        const SAMRAI::hier::BoxOverlap& overlap) const final
        {
            PHARE_LOG_SCOPE(3, "packStream");

            // getDataStreamSize_<true> mean that we want to apply the transformation
            std::size_t expectedSize = getDataStreamSize_<true>(overlap) / sizeof(double);
            std::vector<typename FieldImpl::type> buffer;
            buffer.reserve(expectedSize);

            auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

            SAMRAI::hier::Transformation const& transformation = fieldOverlap.getTransformation();
            if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
            {
                SAMRAI::hier::BoxContainer const& boxContainer
                    = fieldOverlap.getDestinationBoxContainer();
                for (auto const& box : boxContainer)
                {
                    auto const& source = field;
                    SAMRAI::hier::Box sourceBox
                        = Geometry::toFieldBox(getBox(), quantity_, gridLayout);

                    SAMRAI::hier::Box packBox{box};

                    // Since the transformation, allow to transform the source box,
                    // into the destination box space, and that the box in the boxContainer
                    // are in destination space, we have to use the inverseTransform
                    // to get into source space
                    transformation.inverseTransform(packBox);
                    packBox = packBox * sourceBox;

                    internals_.packImpl(buffer, source, packBox, sourceBox);
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
            PHARE_LOG_SCOPE(3, "unpackStream");

            // For unpacking we need to know how much element we will need to
            // extract
            std::size_t expectedSize = getDataStreamSize(overlap) / sizeof(double);

            std::vector<double> buffer;
            buffer.resize(expectedSize, 0.);

            auto& fieldOverlap = dynamic_cast<FieldOverlap const&>(overlap);

            // We flush a portion of the stream on the buffer.
            stream.unpack(buffer.data(), expectedSize);

            SAMRAI::hier::Transformation const& transformation = fieldOverlap.getTransformation();
            if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
            {
                // Here the seek counter will be used to index buffer
                std::size_t seek = 0;

                SAMRAI::hier::BoxContainer const& boxContainer
                    = fieldOverlap.getDestinationBoxContainer();
                for (auto const& box : boxContainer)
                {
                    // For unpackStream, there is no transformation needed, since all the box
                    // are on the destination space

                    auto& source = field;
                    SAMRAI::hier::Box destination
                        = Geometry::toFieldBox(getBox(), quantity_, gridLayout);


                    SAMRAI::hier::Box packBox{box * destination};


                    internals_.unpackImpl(seek, buffer, source, packBox, destination);
                }
            }
        }



        FieldImpl* getPointer() { return &field; }


        static GridLayoutT const& getLayout(SAMRAI::hier::Patch const& patch, int id)
        {
            auto const& patchData = std::dynamic_pointer_cast<FieldData<GridLayoutT, FieldImpl>>(
                patch.getPatchData(id));
            if (!patchData)
            {
                throw std::runtime_error("cannot cast to FieldData");
            }
            return patchData->gridLayout;
        }


        static FieldImpl& getField(SAMRAI::hier::Patch const& patch, int id)
        {
            auto const& patchData = std::dynamic_pointer_cast<FieldData<GridLayoutT, FieldImpl>>(
                patch.getPatchData(id));
            if (!patchData)
            {
                throw std::runtime_error("cannot cast to FieldData");
            }
            return patchData->field;
        }




        GridLayoutT gridLayout;
        FieldImpl field;

    private:
        PhysicalQuantity quantity_; ///! PhysicalQuantity used for this field data




        /*** \brief copy data from the intersection box
         *
         */
        void copy_(SAMRAI::hier::Box const& intersectBox, SAMRAI::hier::Box const& sourceBox,
                   SAMRAI::hier::Box const& destinationBox,
                   [[maybe_unused]] FieldData const& source, FieldImpl const& fieldSource,
                   FieldImpl& fieldDestination)
        {
            // First we represent the intersection that is defined in AMR space to the local space
            // of the source

            SAMRAI::hier::Box localSourceBox{AMRToLocal(intersectBox, sourceBox)};

            // Then we represent the intersection into the local space of the destination
            SAMRAI::hier::Box localDestinationBox{AMRToLocal(intersectBox, destinationBox)};


            // We can finally perform the copy of the element in the correct range
            internals_.copyImpl(localSourceBox, fieldSource, localDestinationBox, fieldDestination);
        }




        void copy_(FieldData const& source, FieldOverlap const& overlap)
        {
            // Here the first step is to get the transformation from the overlap
            // we transform the box from the source, and from the destination
            // from AMR index to FieldData indexes (ie whether or not the quantity is primal
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
                            = Geometry::toFieldBox(source.getBox(), quantity_, source.gridLayout);


                        SAMRAI::hier::Box destinationBox
                            = Geometry::toFieldBox(this->getBox(), quantity_, this->gridLayout);


                        SAMRAI::hier::Box transformedSource{sourceBox};
                        transformation.transform(transformedSource);


                        SAMRAI::hier::Box intersectionBox{box * transformedSource * destinationBox};


                        if (!intersectionBox.empty())
                        {
                            FieldImpl const& sourceField = source.field;
                            FieldImpl& destinationField  = field;

                            copy_(intersectionBox, transformedSource, destinationBox, source,
                                  sourceField, destinationField);
                        }
                    }
                }
            }
            else
            {
                throw std::runtime_error("copy with rotate not implemented");
            }
        }


        /*** \brief Compute the maximum amount of memory needed to hold FieldData information on
         * the specified overlap, this version work on the source, or the destination
         * depending on withTransform parameter
         */
        template<bool withTransform>
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

            // TODO: see FieldDataFactory todo of the same function

            SAMRAI::hier::BoxContainer const& boxContainer
                = fieldOverlap.getDestinationBoxContainer();

            return boxContainer.getTotalSizeOfBoxes() * sizeof(typename FieldImpl::type);
        }


        FieldDataInternals<GridLayoutT, dimension, FieldImpl, PhysicalQuantity> internals_;
    }; // namespace PHARE




    // 1D internals implementation
    template<typename GridLayoutT, typename FieldImpl, typename PhysicalQuantity>
    class FieldDataInternals<GridLayoutT, 1, FieldImpl, PhysicalQuantity>
    {
    public:
        void copyImpl(SAMRAI::hier::Box const& localSourceBox, FieldImpl const& source,
                      SAMRAI::hier::Box const& localDestinationBox, FieldImpl& destination) const
        {
            std::uint32_t xSourceStart = static_cast<std::uint32_t>(localSourceBox.lower(0));
            std::uint32_t xDestinationStart
                = static_cast<std::uint32_t>(localDestinationBox.lower(0));

            std::uint32_t xSourceEnd = static_cast<std::uint32_t>(localSourceBox.upper(0));
            std::uint32_t xDestinationEnd
                = static_cast<std::uint32_t>(localDestinationBox.upper(0));

            for (std::uint32_t xSource = xSourceStart, xDestination = xDestinationStart;
                 xSource <= xSourceEnd && xDestination <= xDestinationEnd;
                 ++xSource, ++xDestination)
            {
                destination(xDestination) = source(xSource);
            }
        }



        void packImpl(std::vector<double>& buffer, FieldImpl const& source,
                      SAMRAI::hier::Box const& overlap, SAMRAI::hier::Box const& sourceBox) const
        {
            int xStart = overlap.lower(0) - sourceBox.lower(0);
            int xEnd   = overlap.upper(0) - sourceBox.lower(0);

            for (int xi = xStart; xi <= xEnd; ++xi)
            {
                buffer.push_back(source(xi));
            }
        }



        void unpackImpl(std::size_t& seek, std::vector<double> const& buffer, FieldImpl& source,
                        SAMRAI::hier::Box const& overlap,
                        SAMRAI::hier::Box const& destination) const
        {
            int xStart = overlap.lower(0) - destination.lower(0);
            int xEnd   = overlap.upper(0) - destination.lower(0);

            for (int xi = xStart; xi <= xEnd; ++xi)
            {
                source(xi) = buffer[seek];
                ++seek;
            }
        }
    };



    // 2D internals implementation
    template<typename GridLayoutT, typename FieldImpl, typename PhysicalQuantity>
    class FieldDataInternals<GridLayoutT, 2, FieldImpl, PhysicalQuantity>
    {
    public:
        void copyImpl(SAMRAI::hier::Box const& localSourceBox, FieldImpl const& source,
                      SAMRAI::hier::Box const& localDestinationBox, FieldImpl& destination) const
        {
            std::uint32_t xSourceStart = static_cast<std::uint32_t>(localSourceBox.lower(0));
            std::uint32_t xDestinationStart
                = static_cast<std::uint32_t>(localDestinationBox.lower(0));

            std::uint32_t xSourceEnd = static_cast<std::uint32_t>(localSourceBox.upper(0));
            std::uint32_t xDestinationEnd
                = static_cast<std::uint32_t>(localDestinationBox.upper(0));

            std::uint32_t ySourceStart = static_cast<std::uint32_t>(localSourceBox.lower(1));
            std::uint32_t yDestinationStart
                = static_cast<std::uint32_t>(localDestinationBox.lower(1));

            std::uint32_t ySourceEnd = static_cast<std::uint32_t>(localSourceBox.upper(1));
            std::uint32_t yDestinationEnd
                = static_cast<std::uint32_t>(localDestinationBox.upper(1));

            for (std::uint32_t xSource = xSourceStart, xDestination = xDestinationStart;
                 xSource <= xSourceEnd && xDestination <= xDestinationEnd;
                 ++xSource, ++xDestination)
            {
                for (std::uint32_t ySource = ySourceStart, yDestination = yDestinationStart;
                     ySource <= ySourceEnd && yDestination <= yDestinationEnd;
                     ++ySource, ++yDestination)
                {
                    destination(xDestination, yDestination) = source(xSource, ySource);
                }
            }
        }




        void packImpl(std::vector<double>& buffer, FieldImpl const& source,
                      SAMRAI::hier::Box const& overlap, SAMRAI::hier::Box const& destination) const

        {
            int xStart = overlap.lower(0) - destination.lower(0);
            int xEnd   = overlap.upper(0) - destination.lower(0);

            int yStart = overlap.lower(1) - destination.lower(1);
            int yEnd   = overlap.upper(1) - destination.lower(1);

            for (int xi = xStart; xi <= xEnd; ++xi)
            {
                for (int yi = yStart; yi <= yEnd; ++yi)
                {
                    buffer.push_back(source(xi, yi));
                }
            }
        }




        void unpackImpl(std::size_t& seek, std::vector<double> const& buffer, FieldImpl& source,
                        SAMRAI::hier::Box const& overlap,
                        SAMRAI::hier::Box const& destination) const
        {
            int xStart = overlap.lower(0) - destination.lower(0);
            int xEnd   = overlap.upper(0) - destination.lower(0);

            int yStart = overlap.lower(1) - destination.lower(1);
            int yEnd   = overlap.upper(1) - destination.lower(1);

            for (int xi = xStart; xi <= xEnd; ++xi)
            {
                for (int yi = yStart; yi <= yEnd; ++yi)
                {
                    source(xi, yi) = buffer[seek];
                    ++seek;
                }
            }
        }
    };



    // 3D internals implementation
    template<typename GridLayoutT, typename FieldImpl, typename PhysicalQuantity>
    class FieldDataInternals<GridLayoutT, 3, FieldImpl, PhysicalQuantity>
    {
    public:
        void copyImpl(SAMRAI::hier::Box const& localSourceBox, FieldImpl const& source,
                      SAMRAI::hier::Box const& localDestinationBox, FieldImpl& destination) const
        {
            std::uint32_t xSourceStart = static_cast<std::uint32_t>(localSourceBox.lower(0));
            std::uint32_t xDestinationStart
                = static_cast<std::uint32_t>(localDestinationBox.lower(0));

            std::uint32_t xSourceEnd = static_cast<std::uint32_t>(localSourceBox.upper(0));
            std::uint32_t xDestinationEnd
                = static_cast<std::uint32_t>(localDestinationBox.upper(0));

            std::uint32_t ySourceStart = static_cast<std::uint32_t>(localSourceBox.lower(1));
            std::uint32_t yDestinationStart
                = static_cast<std::uint32_t>(localDestinationBox.lower(1));

            std::uint32_t ySourceEnd = static_cast<std::uint32_t>(localSourceBox.upper(1));
            std::uint32_t yDestinationEnd
                = static_cast<std::uint32_t>(localDestinationBox.upper(1));

            std::uint32_t zSourceStart = static_cast<std::uint32_t>(localSourceBox.lower(2));
            std::uint32_t zDestinationStart
                = static_cast<std::uint32_t>(localDestinationBox.lower(2));

            std::uint32_t zSourceEnd = static_cast<std::uint32_t>(localSourceBox.upper(2));
            std::uint32_t zDestinationEnd
                = static_cast<std::uint32_t>(localDestinationBox.upper(2));

            for (std::uint32_t xSource = xSourceStart, xDestination = xDestinationStart;
                 xSource <= xSourceEnd && xDestination <= xDestinationEnd;
                 ++xSource, ++xDestination)
            {
                for (std::uint32_t ySource = ySourceStart, yDestination = yDestinationStart;
                     ySource <= ySourceEnd && yDestination <= yDestinationEnd;
                     ++ySource, ++yDestination)
                {
                    for (std::uint32_t zSource = zSourceStart, zDestination = zDestinationStart;
                         zSource <= zSourceEnd && zDestination <= zDestinationEnd;
                         ++zSource, ++zDestination)
                    {
                        destination(xDestination, yDestination, zDestination)
                            = source(xSource, ySource, zSource);
                    }
                }
            }
        }




        void packImpl(std::vector<double>& buffer, FieldImpl const& source,
                      SAMRAI::hier::Box const& overlap, SAMRAI::hier::Box const& destination) const
        {
            int xStart = overlap.lower(0) - destination.lower(0);
            int xEnd   = overlap.upper(0) - destination.lower(0);

            int yStart = overlap.lower(1) - destination.lower(1);
            int yEnd   = overlap.upper(1) - destination.lower(1);

            int zStart = overlap.lower(2) - destination.lower(2);
            int zEnd   = overlap.upper(2) - destination.lower(2);

            for (int xi = xStart; xi <= xEnd; ++xi)
            {
                for (int yi = yStart; yi <= yEnd; ++yi)
                {
                    for (int zi = zStart; zi <= zEnd; ++zi)
                    {
                        buffer.push_back(source(xi, yi, zi));
                    }
                }
            }
        }




        void unpackImpl(std::size_t& seek, std::vector<double> const& buffer, FieldImpl& source,
                        SAMRAI::hier::Box const& overlap,
                        SAMRAI::hier::Box const& destination) const
        {
            int xStart = overlap.lower(0) - destination.lower(0);
            int xEnd   = overlap.upper(0) - destination.lower(0);

            int yStart = overlap.lower(1) - destination.lower(1);
            int yEnd   = overlap.upper(1) - destination.lower(1);

            int zStart = overlap.lower(2) - destination.lower(2);
            int zEnd   = overlap.upper(2) - destination.lower(2);

            for (int xi = xStart; xi <= xEnd; ++xi)
            {
                for (int yi = yStart; yi <= yEnd; ++yi)
                {
                    for (int zi = zStart; zi <= zEnd; ++zi)
                    {
                        source(xi, yi, zi) = buffer[seek];
                        ++seek;
                    }
                }
            }
        }
    };


} // namespace amr
} // namespace PHARE


#endif
