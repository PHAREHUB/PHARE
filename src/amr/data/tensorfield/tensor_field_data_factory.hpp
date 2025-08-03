#ifndef PHARE_SRC_AMR_TENSORFIELD_TENSORFIELD_DATA_FACTORY_HPP
#define PHARE_SRC_AMR_TENSORFIELD_TENSORFIELD_DATA_FACTORY_HPP


#include "core/def/phare_mpi.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include <amr/data/tensorfield/tensor_field_data.hpp>
#include <amr/data/tensorfield/tensor_field_geometry.hpp>

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDataFactory.h>
#include <SAMRAI/tbox/MemoryUtilities.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>

#include <utility>


namespace PHARE::amr
{
template<std::size_t rank, typename GridLayoutT, typename Grid_t, typename PhysicalQuantity>
class TensorFieldDataFactory : public SAMRAI::hier::PatchDataFactory
{
    static constexpr std::size_t n_ghosts
        = GridLayoutT::template nbrGhosts<core::QtyCentering, core::QtyCentering::dual>();

    using tensor_t = typename PhysicalQuantity::template TensorType<rank>;

public:
    static constexpr std::size_t dimension    = GridLayoutT::dimension;
    static constexpr std::size_t interp_order = GridLayoutT::interp_order;


    TensorFieldDataFactory(bool fineBoundaryRepresentsVariable, bool dataLivesOnPatchBorder,
                           std::string const& name, tensor_t qty)
        : SAMRAI::hier::PatchDataFactory(
              SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension(dimension), n_ghosts})
        , fineBoundaryRepresentsVariable_{fineBoundaryRepresentsVariable}
        , dataLivesOnPatchBorder_{dataLivesOnPatchBorder}
        , quantity_{qty}
        , name_{name}
    {
    }




    /*** \brief Clone the current TensorFieldDataFactory
     */
    std::shared_ptr<SAMRAI::hier::PatchDataFactory>
    cloneFactory(SAMRAI::hier::IntVector const& /*ghost*/) final
    {
        return std::make_shared<TensorFieldDataFactory>(fineBoundaryRepresentsVariable_,
                                                        dataLivesOnPatchBorder_, name_, quantity_);
    }




    /*** \brief Given a patch, allocate a TensorFieldData
     * it is expected that this routines will create a functional fieldData
     * (ie with a gridlayout and a Grid_t)
     */
    std ::shared_ptr<SAMRAI::hier::PatchData> allocate(SAMRAI::hier::Patch const& patch) const final
    {
        auto const& domain = patch.getBox();
        SAMRAI::tbox::Dimension dim{dimension};



        // We finally make the TensorFieldData with the correct parameter

        return std::make_shared<TensorFieldData<rank, GridLayoutT, Grid_t, PhysicalQuantity>>(
            domain, SAMRAI::hier::IntVector{dim, n_ghosts}, name_,
            layoutFromPatch<GridLayoutT>(patch), quantity_);
    }




    std::shared_ptr<SAMRAI::hier::BoxGeometry>
    getBoxGeometry(SAMRAI::hier::Box const& box) const final
    {
        // Note : when we create a TensorFieldGeometry, we don't need to have the correct
        // dxdydz, nor the physical origin. All we have to know is the numberCells
        // for the gridlayout, and also we give the box to the TensorFieldGeometry, so that
        // it can use it to get the final box representation.

        std::array<double, dimension> dl;
        std::array<std::uint32_t, dimension> nbCell;
        core::Point<double, dimension> origin;

        for (std::size_t iDim = 0; iDim < dimension; ++iDim)
        {
            dl[iDim]     = 0.01;
            nbCell[iDim] = box.numberCells(iDim);
            origin[iDim] = 0;
        }


        // dumb dl and origin, only nbCell is usefull
        // but TensorFieldGeometry needs to use a gridlayout instance with proper nbrCells
        GridLayoutT gridLayout(dl, nbCell, origin);

        return std::make_shared<TensorFieldGeometry<rank, GridLayoutT, PhysicalQuantity>>(
            box, std::move(gridLayout), quantity_);
    }




    std::size_t getSizeOfMemory(SAMRAI::hier::Box const& box) const final { return 1; }



    bool fineBoundaryRepresentsVariable() const final { return fineBoundaryRepresentsVariable_; }



    bool dataLivesOnPatchBorder() const final { return dataLivesOnPatchBorder_; }



    bool validCopyTo(std::shared_ptr<SAMRAI::hier::PatchDataFactory> const&
                         destinationPatchDataFactory) const final
    {
        auto fieldDataFactory
            = std::dynamic_pointer_cast<TensorFieldDataFactory>(destinationPatchDataFactory);
        return (fieldDataFactory != nullptr);
    }



private:
    bool const fineBoundaryRepresentsVariable_ = false;
    bool const dataLivesOnPatchBorder_         = false;
    tensor_t const quantity_;
    std::string name_;
};


} // namespace PHARE::amr


#endif
