#ifndef PHARE_SRC_AMR_FIELD_FIELD_DATA_FACTORY_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_DATA_FACTORY_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/tbox/MemoryUtilities.h>
#include <SAMRAI/hier/PatchDataFactory.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>

#include "field_data.hpp"

#include <utility>

namespace PHARE
{
namespace amr
{
    template<typename GridLayoutT, typename FieldImpl,
             typename PhysicalQuantity = decltype(std::declval<FieldImpl>().physicalQuantity())>
    /**
     * @brief The FieldDataFactory class
     */
    class FieldDataFactory : public SAMRAI::hier::PatchDataFactory
    {
        static constexpr std::size_t n_ghosts
            = GridLayoutT::template nbrGhosts<core::QtyCentering, core::QtyCentering::dual>();

    public:
        static constexpr std::size_t dimension    = GridLayoutT::dimension;
        static constexpr std::size_t interp_order = GridLayoutT::interp_order;


        FieldDataFactory(bool fineBoundaryRepresentsVariable, bool dataLivesOnPatchBorder,
                         std::string const& name, PhysicalQuantity qty)
            : SAMRAI::hier::PatchDataFactory(
                  SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension(dimension), n_ghosts})
            , fineBoundaryRepresentsVariable_{fineBoundaryRepresentsVariable}
            , dataLivesOnPatchBorder_{dataLivesOnPatchBorder}
            , quantity_{qty}
            , name_{name}
        {
        }




        /*** \brief Clone the current FieldDataFactory
         */
        std::shared_ptr<SAMRAI::hier::PatchDataFactory>
        cloneFactory(SAMRAI::hier::IntVector const& /*ghost*/) final
        {
            return std::make_shared<FieldDataFactory>(fineBoundaryRepresentsVariable_,
                                                      dataLivesOnPatchBorder_, name_, quantity_);
        }




        /*** \brief Given a patch, allocate a FieldData
         * it is expected that this routines will create a functional fieldData
         * (ie with a gridlayout and a FieldImpl)
         */
        std ::shared_ptr<SAMRAI::hier::PatchData>
        allocate(SAMRAI::hier::Patch const& patch) const final
        {
            auto const& domain = patch.getBox();
            SAMRAI::tbox::Dimension dim{dimension};



            // We finally make the FieldData with the correct parameter

            return std::make_shared<FieldData<GridLayoutT, FieldImpl>>(
                domain, SAMRAI::hier::IntVector{dim, n_ghosts}, name_,
                layoutFromPatch<GridLayoutT>(patch), quantity_);
        }




        std::shared_ptr<SAMRAI::hier::BoxGeometry>
        getBoxGeometry(SAMRAI::hier::Box const& box) const final
        {
            // Note : when we create a FieldGeometry, we don't need to have the correct
            // dxdydz, nor the physical origin. All we have to know is the numberCells
            // for the gridlayout, and also we give the box to the FieldGeometry, so that
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
            // but FieldGeometry needs to use a gridlayout instance with proper nbrCells
            GridLayoutT gridLayout(dl, nbCell, origin);

            return std::make_shared<FieldGeometry<GridLayoutT, PhysicalQuantity>>(
                box, std::move(gridLayout), quantity_);
        }




        std::size_t getSizeOfMemory(SAMRAI::hier::Box const& box) const final
        {
            // TODO: this calculus assumes that we don't need more memory than
            //       alignedMemory(nx*ny*nz*sizeof(double)) + alignedMemory(baseSize)

            std::array<double, dimension> dl;
            std::array<std::uint32_t, dimension> nbCell;
            core::Point<double, dimension> origin;

            for (std::size_t iDim = 0; iDim < dimension; ++iDim)
            {
                dl[iDim]     = 0.01; // some value that is not used anyway
                origin[iDim] = 0;
                nbCell[iDim] = box.numberCells(iDim);
            }

            std::size_t const baseField
                = SAMRAI::tbox::MemoryUtilities::align(sizeof(FieldData<GridLayoutT, FieldImpl>));

            GridLayoutT gridLayout{dl, nbCell, origin};


            auto const& allocSize = gridLayout.allocSize(quantity_);

            std::size_t data = 1;
            for (auto nCell : allocSize)
            {
                data *= nCell;
            }

            data *= sizeof(typename FieldImpl::type);



            return baseField + SAMRAI::tbox::MemoryUtilities::align(data);
        }




        bool fineBoundaryRepresentsVariable() const final
        {
            return fineBoundaryRepresentsVariable_;
        }




        bool dataLivesOnPatchBorder() const final { return dataLivesOnPatchBorder_; }



        bool validCopyTo(std::shared_ptr<SAMRAI::hier::PatchDataFactory> const&
                             destinationPatchDataFactory) const final
        {
            auto fieldDataFactory
                = std::dynamic_pointer_cast<FieldDataFactory>(destinationPatchDataFactory);
            return (fieldDataFactory != nullptr);
        }



    private:
        bool const fineBoundaryRepresentsVariable_ = false;
        bool const dataLivesOnPatchBorder_         = false;
        PhysicalQuantity const quantity_;
        std::string name_;
    };


} // namespace amr

} // namespace PHARE


#endif
