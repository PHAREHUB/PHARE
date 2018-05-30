#ifndef PHARE_SRC_AMR_FIELD_FIELD_DATA_H
#define PHARE_SRC_AMR_FIELD_FIELD_DATA_H

#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/tbox/MemoryUtilities.h>
#include <utility>

#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"

#include "field_geometry.h"

namespace PHARE
{
// We use another class here so that we can specialize specifics function: copy , pack , unpack
// on the dimension and we don't want to loose non specialized function related to SAMRAI interface
template<typename GridLayoutImpl, std::size_t dim, typename FieldImpl,
         typename PhysicalQuantity = decltype(std::declval<FieldImpl>().physicalQuantity())>
class FieldDataInternals
{
};



template<typename GridLayoutImpl, typename FieldImpl,
         typename PhysicalQuantity = decltype(std::declval<FieldImpl>().physicalQuantity())>
class FieldData : public SAMRAI::hier::PatchData
{
public:
    static constexpr std::size_t dimension    = GridLayoutImpl::dimension;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;

    /*** \brief Construct a FieldData from information associated to a patch
     *
     * It will create a GridLayout from parameters given by FieldDataFactory
     * From the freshly created GridLayout, it will create a Field with the correct
     * number of cells in each needed directions
     */
    FieldData(SAMRAI::hier::Box const& domain, SAMRAI::hier::IntVector const& ghost,
              std::string name, std::array<double, dimension> const& dl,
              std::array<uint32, dimension> const& nbrCells, Point<double, dimension> const& origin,
              PhysicalQuantity qty)
        : SAMRAI::hier::PatchData(domain, ghost)
        , gridLayout{dl, nbrCells, origin}
        , field(name, qty, gridLayout.allocSize(qty))
        , quantity_{qty}
    {
    }




    FieldData()                 = delete;
    FieldData(FieldData const&) = delete;
    FieldData(FieldData&&)      = default;

    FieldData& operator=(FieldData const&) = delete;



    /*** \brief Copy information from another FieldData where data overlap
     *
     *    The data will be copied from the interior and ghost of the source to the interior and
     *    ghost of the destination, where there is an overlap in the underlying index space
     */
    void copy(const SAMRAI::hier::PatchData& source) final
    {
        // After checking that source and *this have the same number of dimension
        // We will try to cast source as a FieldData, if it succeed we can continue
        // and perform the copy. Otherwise we call copy2 that will simply throw a runtime
        // error

        TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source);

        auto fieldSource = dynamic_cast<FieldData const*>(&source);


        if (fieldSource != nullptr)
        {
            // First step is to translate the AMR box into GridLayout index space
            // to accomplish that we get the interior box, from the FieldData.
            // and we call toFieldBox with the parameter withGhost = true.
            // note that we could have stored the ghost box of the field data at
            // creation

            SAMRAI::hier::Box sourceBox
                = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
                    fieldSource->getBox(), quantity_, fieldSource->gridLayout);


            SAMRAI::hier::Box destinationBox
                = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
                    this->getBox(), quantity_, this->gridLayout);

            // Given the two box in correct space we just have to intersect them
            SAMRAI::hier::Box intersectionBox = sourceBox * destinationBox;

            if (!intersectionBox.empty())
            {
                auto const& sourceField = fieldSource->field;
                auto& destinationField  = field;

                // We can copy field from the source to the destination on the correct region
                copy_(intersectionBox, sourceBox, destinationBox, *fieldSource, sourceField,
                      destinationField);
            }
        }
        else
        {
            // If we go there, then the PatchData is not a FieldData
            // We cannot be in this case since FieldData is final
            // hence copy2 will directly throw
            // Note that we still have to implement copy2 and to call it, since it is
            // what is expected from SAMRAI
            source.copy2(*this);
        }
    }




    /*** \brief This form should not be called since we cannot derive from FieldData
     * since FieldData is a final implementation of PatchData
     */
    void copy2(SAMRAI::hier::PatchData& destination) const final
    {
        throw std::runtime_error("Error cannot cast the PatchData to FieldData");
    }




    /*** \brief Copy data from the source into the destination using the designated overlap
     * descriptor.
     *
     *   The overlap will contain information on which part have to be copied, and where it will be
     * fill, in the destination. It will also give the necessary transformation to apply to the
     * source, to perform the copy (ie : translation for periodics condition)
     */
    void copy(const SAMRAI::hier::PatchData& source, const SAMRAI::hier::BoxOverlap& overlap) final
    {
        auto fieldSource  = dynamic_cast<FieldData const*>(&source);
        auto fieldOverlap = dynamic_cast<FieldOverlap<dimension> const*>(&overlap);


        // So here we just check that the PatchData is a FieldData, and that this is a correct
        // FieldOverlap
        if ((fieldSource != nullptr) && (fieldOverlap != nullptr))
        {
            // If it is the case,we delegate the copy on another function
            // that take directly the FieldData, and FieldOverlap

            copy_(*fieldSource, *fieldOverlap);
        }
        else
        {
            source.copy2(*this, overlap);
        }
    }




    /*** \brief This form should not be called since we cannot derive from FieldData
     */
    void copy2(SAMRAI::hier::PatchData& destination,
               const SAMRAI::hier::BoxOverlap& overlap) const final
    {
        throw std::runtime_error("Error cannot cast the PatchData to FieldData");
    }




    /*** \brief Determines whether the patch data subclass can estimate the necessary stream size
     * using only index space information.
     *
     * The return value is true since that for a corresponding domain, there is a fixed
     * number of elements in the field depending on the PhysicalQuantity and the Layout used
     */
    bool canEstimateStreamSizeFromBox() const final { return true; }



    /*** \brief Compute the maximum amount of memory needed to hold FieldData information on
     * the specified overlap
     */
    size_t getDataStreamSize(const SAMRAI::hier::BoxOverlap& overlap) const final
    {
        // The idea here is to tell SAMRAI the maximum memory will be used by our type
        // on a given region.
        // this version assume that there is no transformation needed to be applied.
        // It will be the case where it is used with unpackStream.

        FieldOverlap<dimension> const* fieldOverlap
            = dynamic_cast<FieldOverlap<dimension> const*>(&overlap);
        TBOX_ASSERT(fieldOverlap != nullptr);

        size_t totalSize = 0;

        if (fieldOverlap->isOverlapEmpty())
        {
            return 0u;
        }

        // TODO: see FieldDataFactory todo of the same function

        SAMRAI::hier::BoxContainer const& boxContainer = fieldOverlap->getDestinationBoxContainer();

        for (auto const& box : boxContainer)
        {
            // We compute the intersection between the box contained in the overlap
            // with the ghostBox of the fieldData (toFieldBox , withGhost=true default parameter)
            SAMRAI::hier::Box finalBox{box};
            finalBox = finalBox
                       * FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
                             getBox(), quantity_, gridLayout);

            size_t size = 1;

            for (uint32 iDir = 0; iDir < dimension; ++iDir)
            {
                size *= finalBox.numberCells(iDir);
            }

            // size *= sizeof(double);
            // At the end we will make sure that the size correspond to an aligned memory
            // since it will be the max memory used
            totalSize += size;
        }
        totalSize
            = SAMRAI::tbox::MemoryUtilities::align(totalSize * sizeof(typename FieldImpl::type));
        return totalSize;
    }




    /*** \brief Serialize the data contained in the field data on the region covered by the overlap,
     * and put it on the stream.
     */
    void packStream(SAMRAI::tbox::MessageStream& stream,
                    const SAMRAI::hier::BoxOverlap& overlap) const final
    {
        // TODO: getDataStreamSize can only get the size of data when unpacking
        // so we cannot use it with packStream.
        // Note that we should refactor into two func, one with the overlap transformation
        // and the other without, so that getDataStreamSize will call the one without the
        // transformation and here we will call the one with the transformation

        /* TODO size_t expectedSize = getDataStreamSize(overlap) / sizeof(double); */
        std::vector<typename FieldImpl::type> buffer;
        /* buffer.reserve(expectedSize); */

        auto fieldOverlap = dynamic_cast<FieldOverlap<dimension> const*>(&overlap);
        TBOX_ASSERT(fieldOverlap != nullptr);

        SAMRAI::hier::Transformation const& transformation = fieldOverlap->getTransformation();
        if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
        {
            SAMRAI::hier::BoxContainer const& boxContainer
                = fieldOverlap->getDestinationBoxContainer();
            for (auto const& box : boxContainer)
            {
                auto const& source = field;
                SAMRAI::hier::Box sourceBox
                    = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
                        getBox(), quantity_, gridLayout);

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




    /*** \brief Unserialize data contained on the stream, that comes from a region covered by the
     * overlap, and fill the data where is needed.
     */
    void unpackStream(SAMRAI::tbox::MessageStream& stream,
                      const SAMRAI::hier::BoxOverlap& overlap) final
    {
        // For unpacking we need to know how much element we will need to
        // extract
        size_t expectedSize = getDataStreamSize(overlap) / sizeof(double);
        // Here the seek counter will be used to index buffer
        size_t seek = 0;
        std::vector<double> buffer;
        buffer.resize(expectedSize, 0.);

        auto fieldOverlap = dynamic_cast<FieldOverlap<dimension> const*>(&overlap);
        TBOX_ASSERT(fieldOverlap != nullptr);


        // We flush a portion of the stream on the buffer.
        stream.unpack(buffer.data(), expectedSize);

        SAMRAI::hier::Transformation const& transformation = fieldOverlap->getTransformation();
        if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
        {
            SAMRAI::hier::BoxContainer const& boxContainer
                = fieldOverlap->getDestinationBoxContainer();
            for (auto const& box : boxContainer)
            {
                // For unpackStream, there is no transformation needed, since all the box
                // are on the destination space

                auto& source = field;
                SAMRAI::hier::Box destination
                    = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
                        getBox(), quantity_, gridLayout);


                SAMRAI::hier::Box packBox{box * destination};


                internals_.unpackImpl(seek, buffer, source, packBox, destination);
            }
        }
    }




    GridLayout<GridLayoutImpl> gridLayout;
    FieldImpl field;

private:
    PhysicalQuantity quantity_; ///! PhysicalQuantity used for this field data




    /*** \brief copy data from the intersection box
     *
     */
    void copy_(SAMRAI::hier::Box const& intersectBox, SAMRAI::hier::Box const& sourceBox,
               SAMRAI::hier::Box const& destinationBox, FieldData const& source,
               FieldImpl const& fieldSource, FieldImpl& fieldDestination)
    {
        // First we represent the intersection that is defined in AMR space to the local space of
        // the source
        SAMRAI::hier::Box localSourceBox(intersectBox.getDim());
        localSourceBox.setLower(intersectBox.lower() - sourceBox.lower());
        localSourceBox.setUpper(intersectBox.upper() - sourceBox.lower());

        // Then we represent the intersection into the local space of the destination
        SAMRAI::hier::Box localDestinationBox(intersectBox.getDim());
        localDestinationBox.setLower(intersectBox.lower() - destinationBox.lower());
        localDestinationBox.setUpper(intersectBox.upper() - destinationBox.lower());

        // We can finnaly perform the copy of the element in the correct range
        internals_.copyImpl(localSourceBox, fieldSource, localDestinationBox, fieldDestination);
    }




    void copy_(FieldData const& source, FieldOverlap<dimension> const& overlap)
    {
        // Here the first step is to get the transformation from the overlap
        // then we have two case (We  only allow translation transformation):
        // case one: we  have a null translation
        //          |
        //          -> In this case we transform the box from the source, and from the destination
        //             from AMR index to FieldData indexes (ie whether or not the quantity is primal
        //             or not), and we also consider the ghost. After that we compute the
        //             intersection with the source box, the destinationBox, and the box from the
        //             destinationBoxContainer.
        //
        //
        // case two: we have a non null translation
        //          |
        //          -> same as case one, except that we apply the transformation to the source box


        SAMRAI::hier::Transformation const& transformation = overlap.getTransformation();

        if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE)
        {
            SAMRAI::hier::BoxContainer const& boxList = overlap.getDestinationBoxContainer();

            SAMRAI::hier::IntVector const zeroOffset{
                SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension{dimension})};

            if (transformation.getOffset() == zeroOffset
                && transformation.getBeginBlock() == transformation.getEndBlock())
            {
                for (auto const& box : boxList)
                {
                    SAMRAI::hier::Box sourceBox
                        = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
                            source.getBox(), quantity_, source.gridLayout);
                    SAMRAI::hier::Box destinationBox
                        = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
                            this->getBox(), quantity_, this->gridLayout);

                    SAMRAI::hier::Box intersectionBox{box * sourceBox * destinationBox};

                    if (!intersectionBox.empty())
                    {
                        FieldImpl const& sourceField = source.field;
                        FieldImpl& destinationField  = field;
                        copy_(intersectionBox, sourceBox, destinationBox, source, sourceField,
                              destinationField);
                    }
                }
            }
            else
            {
                for (auto const& box : boxList)
                {
                    SAMRAI::hier::Box sourceBox
                        = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
                            source.getBox(), quantity_, source.gridLayout);
                    SAMRAI::hier::Box destinationBox
                        = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
                            this->getBox(), quantity_, this->gridLayout);
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



    FieldDataInternals<GridLayoutImpl, dimension, FieldImpl, PhysicalQuantity> internals_;
}; // namespace PHARE




// 1D internals implementation
template<typename GridLayoutImpl, typename FieldImpl, typename PhysicalQuantity>
class FieldDataInternals<GridLayoutImpl, 1, FieldImpl, PhysicalQuantity>
{
public:
    void copyImpl(SAMRAI::hier::Box const& localSourceBox, FieldImpl const& source,
                  SAMRAI::hier::Box const& localDestinationBox, FieldImpl& destination) const
    {
        uint32 xSourceStart      = localSourceBox.lower(0);
        uint32 xDestinationStart = localDestinationBox.lower(0);

        uint32 xSourceEnd      = localSourceBox.upper(0);
        uint32 xDestinationEnd = localDestinationBox.upper(0);

        for (uint32 xSource = xSourceStart, xDestination = xDestinationStart;
             xSource <= xSourceEnd && xDestination <= xDestinationEnd; ++xSource, ++xDestination)
        {
            destination(xDestination) = source(xSource);
        }
    }




    void packImpl(std::vector<double>& buffer, FieldImpl const& source,
                  SAMRAI::hier::Box const& overlap, SAMRAI::hier::Box const& destination) const
    {
        int xStart = overlap.lower(0) - destination.lower(0);
        int xEnd   = overlap.upper(0) - destination.lower(0);

        for (int xi = xStart; xi <= xEnd; ++xi)
        {
            buffer.push_back(source(xi));
        }
    }




    void unpackImpl(size_t& seek, std::vector<double> const& buffer, FieldImpl& source,
                    SAMRAI::hier::Box const& overlap, SAMRAI::hier::Box const& destination) const
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
template<typename GridLayoutImpl, typename FieldImpl, typename PhysicalQuantity>
class FieldDataInternals<GridLayoutImpl, 2, FieldImpl, PhysicalQuantity>
{
public:
    void copyImpl(SAMRAI::hier::Box const& localSourceBox, FieldImpl const& source,
                  SAMRAI::hier::Box const& localDestinationBox, FieldImpl& destination) const
    {
        uint32 xSourceStart      = localSourceBox.lower(0);
        uint32 xDestinationStart = localDestinationBox.lower(0);

        uint32 xSourceEnd      = localSourceBox.upper(0);
        uint32 xDestinationEnd = localDestinationBox.upper(0);

        uint32 ySourceStart      = localSourceBox.lower(1);
        uint32 yDestinationStart = localDestinationBox.lower(1);

        uint32 ySourceEnd      = localSourceBox.upper(1);
        uint32 yDestinationEnd = localDestinationBox.upper(1);

        for (uint32 xSource = xSourceStart, xDestination = xDestinationStart;
             xSource <= xSourceEnd && xDestination <= xDestinationEnd; ++xSource, ++xDestination)
        {
            for (uint32 ySource = ySourceStart, yDestination = yDestinationStart;
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




    void unpackImpl(size_t& seek, std::vector<double> const& buffer, FieldImpl& source,
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
                source(xi, yi) = buffer[seek];
                ++seek;
            }
        }
    }
};



// 3D internals implementation
template<typename GridLayoutImpl, typename FieldImpl, typename PhysicalQuantity>
class FieldDataInternals<GridLayoutImpl, 3, FieldImpl, PhysicalQuantity>
{
public:
    void copyImpl(SAMRAI::hier::Box const& localSourceBox, FieldImpl const& source,
                  SAMRAI::hier::Box const& localDestinationBox, FieldImpl& destination) const
    {
        uint32 xSourceStart      = localSourceBox.lower(0);
        uint32 xDestinationStart = localDestinationBox.lower(0);

        uint32 xSourceEnd      = localSourceBox.upper(0);
        uint32 xDestinationEnd = localDestinationBox.upper(0);

        uint32 ySourceStart      = localSourceBox.lower(1);
        uint32 yDestinationStart = localDestinationBox.lower(1);

        uint32 ySourceEnd      = localSourceBox.upper(1);
        uint32 yDestinationEnd = localDestinationBox.upper(1);

        uint32 zSourceStart      = localSourceBox.lower(2);
        uint32 zDestinationStart = localDestinationBox.lower(2);

        uint32 zSourceEnd      = localSourceBox.upper(2);
        uint32 zDestinationEnd = localDestinationBox.upper(2);

        for (uint32 xSource = xSourceStart, xDestination = xDestinationStart;
             xSource <= xSourceEnd && xDestination <= xDestinationEnd; ++xSource, ++xDestination)
        {
            for (uint32 ySource = ySourceStart, yDestination = yDestinationStart;
                 ySource <= ySourceEnd && yDestination <= yDestinationEnd;
                 ++ySource, ++yDestination)
            {
                for (uint32 zSource = zSourceStart, zDestination = zDestinationStart;
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




    void unpackImpl(size_t& seek, std::vector<double> const& buffer, FieldImpl& source,
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
                    source(xi, yi, zi) = buffer[seek];
                    ++seek;
                }
            }
        }
    }
};




} // namespace PHARE


#endif
