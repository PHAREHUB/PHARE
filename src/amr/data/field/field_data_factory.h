#ifndef PHARE_SRC_AMR_FIELD_FIELD_DATA_FACTORY_H
#define PHARE_SRC_AMR_FIELD_FIELD_DATA_FACTORY_H

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDataFactory.h>
#include <SAMRAI/tbox/MemoryUtilities.h>

#include <utility>

#include "data/coarsening/field_coarsen.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "field_data.h"

namespace PHARE
{
template<typename GridLayoutImpl, typename FieldImpl,
         typename PhysicalQuantity = decltype(std::declval<FieldImpl>().physicalQuantity())>
class FieldDataFactory : public SAMRAI::hier::PatchDataFactory
{
public:
    static constexpr std::size_t dimension    = GridLayoutImpl::dimension;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;


    FieldDataFactory(bool fineBoundaryRepresentsVariable, bool dataLivesOnPatchBorder,
                     std::string const& name, PhysicalQuantity qty)
        : SAMRAI::hier::PatchDataFactory(
              SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension(dimension)))
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
        (void)ghost;
        return std::make_shared<FieldDataFactory>(fineBoundaryRepresentsVariable_,
                                                  dataLivesOnPatchBorder_, name_, quantity_);
    }




    /*** \brief Given a patch, allocate a FieldData
     * it is expected that this routines will create a functional fieldData
     * (ie with a gridlayout and a FieldImpl)
     */
    std ::shared_ptr<SAMRAI::hier::PatchData> allocate(SAMRAI::hier::Patch const& patch) const final
    {
        SAMRAI::tbox::Dimension const dim{dimension};
        //  We get geometry information from the patch, such as meshSize, and physical origin
        auto patchGeom = std::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(
            patch.getPatchGeometry());
        Point<double, dimension> origin;

        std::array<double, dimension> dl;

        if (patchGeom != nullptr)
        {
            auto pOrigin = patchGeom->getXLower();
            auto pDl     = patchGeom->getDx();

            for (std::size_t iDim = 0; iDim < dimension; ++iDim)
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

        std::array<uint32, dimension> nbrCell;

        for (std::size_t iDim = 0; iDim < dimension; ++iDim)
        {
            nbrCell[iDim] = static_cast<uint32>(domain.numberCells(iDim));
        }


        // We finnaly make the FieldData with the correct parameter

        return std::make_shared<FieldData<GridLayoutImpl, FieldImpl>>(
            domain, SAMRAI::hier::IntVector::getZero(dim), name_, dl, nbrCell, origin, quantity_);
    }




    std::shared_ptr<SAMRAI::hier::BoxGeometry>
    getBoxGeometry(SAMRAI::hier::Box const& box) const final
    {
        // Note : when we create a FieldGeometry, we don't need to have the correct
        // dxdydz, nor the physical origin. All we have to know is the numberCells
        // for the gridlayout, and also we give the box to the FieldGeometry, so that
        // it can use it to get the final box representation.

        std::array<double, dimension> dl;
        std::array<uint32, dimension> nbCell;
        Point<double, dimension> origin;

        for (std::size_t iDim = 0; iDim < dimension; ++iDim)
        {
            dl[iDim]     = 0.01;
            nbCell[iDim] = box.numberCells(iDim);
            origin[iDim] = 0;
        }


        GridLayout<GridLayoutImpl> gridLayout(dl, nbCell, origin);

        return std::make_shared<FieldGeometry<GridLayoutImpl, PhysicalQuantity>>(
            box, std::move(gridLayout), quantity_);
    }




    std::size_t getSizeOfMemory(SAMRAI::hier::Box const& box) const final
    {
        // TODO: this calculus assumes that we don't need more memory than
        //       alignedMemory(nx*ny*nz*sizeof(double)) + alignedMemory(baseSize)

        std::array<double, dimension> dl;
        std::array<uint32, dimension> nbCell;
        Point<double, dimension> origin;

        for (std::size_t iDim = 0; iDim < dimension; ++iDim)
        {
            dl[iDim]     = 0.01; // some value that is not used anyway
            origin[iDim] = 0;
            nbCell[iDim] = box.numberCells(iDim);
        }

        const std::size_t baseField
            = SAMRAI::tbox::MemoryUtilities::align(sizeof(FieldData<GridLayoutImpl, FieldImpl>));

        GridLayout<GridLayoutImpl> gridLayout{dl, nbCell, origin};


        auto const& allocSize = gridLayout.allocSize(quantity_);

        std::size_t data = 1;
        for (auto nCell : allocSize)
        {
            data *= nCell;
        }

        data *= sizeof(typename FieldImpl::type);



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
    std::string name_;
};




} // namespace PHARE


#endif
