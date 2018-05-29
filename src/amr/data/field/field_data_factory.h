#ifndef PHARE_SRC_AMR_FIELD_FIELD_DATA_FACTORY_H
#define PHARE_SRC_AMR_FIELD_FIELD_DATA_FACTORY_H

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDataFactory.h>
#include <SAMRAI/tbox/MemoryUtilities.h>

#include <utility>

#include "data/grid/gridlayout.h"
#include "field_data.h"

namespace PHARE
{
template<Layout layout, std::size_t dim, std::size_t interpOrder, typename FieldImpl,
         typename PhysicalQuantity = decltype(std::declval<FieldImpl>().physicalQuantity())>
class FieldDataFactory : public SAMRAI::hier::PatchDataFactory
{
public:
    FieldDataFactory(bool fineBoundaryRepresentsVariable, bool dataLivesOnPatchBorder,
                     std::string const& name, PhysicalQuantity qty)
        : SAMRAI::hier::PatchDataFactory(
              SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension(dim)))
        , fineBoundaryRepresentsVariable_{fineBoundaryRepresentsVariable}
        , dataLivesOnPatchBorder_{dataLivesOnPatchBorder}
        , quantity_{qty}
        , name_{name}
    {
    }




    /*** \brief Clone the current FieldDataFactory
     */
    std::shared_ptr<SAMRAI::hier::PatchDataFactory>
    cloneFactory(SAMRAI::hier::IntVector const& ghost) final
    {
        return std::make_shared<FieldDataFactory>(fineBoundaryRepresentsVariable_,
                                                  dataLivesOnPatchBorder_, name_, quantity_);
    }




    /*** \brief Given a patch, allocate a FieldData
     * it is expected that this routines will create a functional fieldData
     * (ie with a gridlayout and a FieldImpl)
     */
    std ::shared_ptr<SAMRAI::hier::PatchData> allocate(SAMRAI::hier::Patch const& patch) const final
    {
        SAMRAI::tbox::Dimension const dimension{dim};
        //  We get geometry information from the patch, such as meshSize, and physical origin
        auto patchGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            patch.getPatchGeometry());
        Point<double, dim> origin;

        std::array<double, dim> dl;

        if (patchGeom != nullptr)
        {
            auto pOrigin = patchGeom->getXLower();
            auto pDl     = patchGeom->getDx();

            for (std::size_t iDim = 0; iDim < dim; ++iDim)
            {
                origin[iDim] = pOrigin[iDim];
                dl[iDim]     = pDl[iDim];
            }
        }
        else
        {
            // in case that the patch does not have a CartesianPatchGeometry
            // the gridlayout will most likely throw at the construction
            // so we may throw here instead
            throw std::runtime_error(
                "The geometry on the patch is not set, please verify your configuration");
        }

        SAMRAI::hier::Box domain = patch.getBox();

        std::array<uint32, dim> nbrCell;

        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            nbrCell[iDim] = static_cast<uint32>(domain.numberCells(iDim));
        }


        // We finnaly make the FieldData with the correct parameter

        return std::make_shared<FieldData<layout, dim, interpOrder, FieldImpl, PhysicalQuantity>>(
            domain, SAMRAI::hier::IntVector::getZero(dimension), name_, dl, nbrCell, layoutName_,
            origin, quantity_);
    }




    std::shared_ptr<SAMRAI::hier::BoxGeometry>
    getBoxGeometry(SAMRAI::hier::Box const& box) const final
    {
        // Note : when we create a FieldGeometry, we don't need to have the correct
        // dxdydz, nor the physical origin. All we have to know is the numberCells
        // for the gridlayout, and also we give the box to the FieldGeometry, so that
        // it can use it to get the final box representation.

        std::array<double, dim> dl;
        std::array<uint32, dim> nbCell;
        Point<double, dim> origin;

        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            dl[iDim]     = 0.01;
            nbCell[iDim] = box.numberCells(iDim);
            origin[iDim] = 0;
        }


        GridLayout<layout, dim> gridLayout(dl, nbCell, name_, origin, interpOrder);

        return std::make_shared<FieldGeometry<dim, decltype(gridLayout), PhysicalQuantity>>(
            box, std::move(gridLayout), quantity_);
    }




    std::size_t getSizeOfMemory(SAMRAI::hier::Box const& box) const final
    {
        // TODO: this calculus assumes that we don't need more memory than
        //       alignedMemory(nx*ny*nz*sizeof(double)) + alignedMemory(baseSize)

        std::array<double, dim> dl;
        std::array<uint32, dim> nbCell;
        Point<double, dim> origin;

        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            dl[iDim]     = 0.01; // some value that is not used anyway
            origin[iDim] = 0;
            nbCell[iDim] = box.numberCells(iDim);
        }

        const std::size_t baseField = SAMRAI::tbox::MemoryUtilities::align(
            sizeof(FieldData<layout, dim, interpOrder, FieldImpl, PhysicalQuantity>));

        GridLayout<layout, dim> gridLayout(dl, nbCell, name_, origin, interpOrder);


        auto const& allocSize = gridLayout.allocSize(quantity_);

        std::size_t data = 1;
        for (auto nCell : allocSize)
        {
            data *= nCell;
        }

        data *= sizeof(double); // TODO ----> FieldImpl::type



        return baseField + SAMRAI::tbox::MemoryUtilities::align(data);
    }




    bool fineBoundaryRepresentsVariable() const final { return fineBoundaryRepresentsVariable_; }




    bool dataLivesOnPatchBorder() const final { return dataLivesOnPatchBorder_; }



    bool validCopyTo(std::shared_ptr<SAMRAI::hier::PatchDataFactory> const&
                         destinationPatchDataFactory) const final
    {
        auto fieldDataFactory
            = std::dynamic_pointer_cast<FieldDataFactory>(destinationPatchDataFactory);
        return (fieldDataFactory != nullptr);
    }



private:
    bool const fineBoundaryRepresentsVariable_;
    bool const dataLivesOnPatchBorder_;
    PhysicalQuantity const quantity_;
    std::string const layoutName_{
        "yee"}; // not sure that gridlayout still needs a name // TODO  indeed budy
    std::string name_;
};




} // namespace PHARE


#endif
