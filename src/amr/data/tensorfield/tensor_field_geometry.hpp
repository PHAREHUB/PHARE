#ifndef PHARE_SRC_AMR_TENSORFIELD_TENSORFIELD_GEOMETRY_HPP
#define PHARE_SRC_AMR_TENSORFIELD_TENSORFIELD_GEOMETRY_HPP


#include "amr/data/field/field_geometry.hpp"
#include "amr/data/tensorfield/tensor_field_overlap.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/utilities/types.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "amr/data/field/field_overlap.hpp"

#include <SAMRAI/hier/Box.h>
#include "SAMRAI/hier/IntVector.h"
#include <SAMRAI/hier/BoxGeometry.h>


#include <array>
#include <cassert>
#include <type_traits>

namespace PHARE::amr
{


template<std::size_t dimension, std::size_t rank = 1>
class TensorFieldGeometryBase : public SAMRAI::hier::BoxGeometry
{
    static constexpr auto N = core::detail::tensor_field_dim_from_rank<rank>();

public:
    virtual ~TensorFieldGeometryBase() {}
    TensorFieldGeometryBase(SAMRAI::hier::Box const& patch_box,
                            SAMRAI::hier::Box const& ghostTensorFieldBox,
                            SAMRAI::hier::Box const& interiorTensorFieldBox,
                            std::array<core::QtyCentering, dimension> const& centerings)
        : patchBox{patch_box}
        , ghostTensorFieldBox_{ghostTensorFieldBox}
        , interiorTensorFieldBox_{interiorTensorFieldBox}
        , centerings_{centerings}
    {
    }

    auto const& interiorTensorFieldBox() const { return interiorTensorFieldBox_; }

    SAMRAI::hier::Box const patchBox;

protected:
    std::array<SAMRAI::hier::Box, N> const ghostTensorFieldBox_;
    std::array<SAMRAI::hier::Box, N> const interiorTensorFieldBox_;
    std::array<core::QtyCentering, dimension> const centerings_;
};

template<std::size_t rank, typename GridLayoutT, typename PhysicalQuantity>
class TensorFieldGeometry : public TensorFieldGeometryBase<GridLayoutT::dimension>
{
    using tensor_t          = typename PhysicalQuantity::template TensorType<rank>;
    static constexpr auto N = core::detail::tensor_field_dim_from_rank<rank>();

public:
    using Super                               = TensorFieldGeometryBase<GridLayoutT::dimension>;
    static constexpr std::size_t dimension    = GridLayoutT::dimension;
    static constexpr std::size_t interp_order = GridLayoutT::interp_order;


    TensorFieldGeometry(SAMRAI::hier::Box const& box, GridLayoutT const& layout, tensor_t const qty)
        : Super(box,
                toFieldBox(SAMRAI::hier::Box::grow(
                               box, SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dimension},
                                                            GridLayoutT::nbrGhosts()}),
                           qty, layout),
                toFieldBox(box, qty, layout),
                ConstArray<core::QtyCentering, dimension>(core::QtyCentering::primal))
        , layout_{layout}
        , quantity_{qty}
    {
    }




    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& destinationGeometry,
                     SAMRAI::hier::BoxGeometry const& sourceGeometry,
                     SAMRAI::hier::Box const& sourceMask, SAMRAI::hier::Box const& fillBox,
                     bool const overwriteInterior, SAMRAI::hier::Transformation const& sourceOffset,
                     [[maybe_unused]] bool const retry,
                     SAMRAI::hier::BoxContainer const& destinationRestrictBoxes
                     = SAMRAI::hier::BoxContainer{}) const final
    {
        auto& destinationCast = dynamic_cast<TensorFieldGeometry const&>(destinationGeometry);
        auto& sourceCast      = dynamic_cast<TensorFieldGeometry const&>(sourceGeometry);
        return doOverlap_(destinationCast, sourceCast, sourceMask, fillBox, overwriteInterior,
                          sourceOffset, destinationRestrictBoxes);
    }




    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    setUpOverlap(SAMRAI::hier::BoxContainer const& boxes,
                 SAMRAI::hier::Transformation const& offset) const final
    {
        auto qts = PhysicalQuantity::componentsQuantities(quantity_);

        auto overlaps = core::for_N<N, core::for_N_R_mode::make_array>([&](auto i) {
            SAMRAI::hier::BoxContainer destinationBoxes;

            auto const qty = qts[i];

            for (auto& box : boxes)
            {
                core::GridLayout const layout = layoutFromBox(box, layout_);
                destinationBoxes.push_back(
                    FieldGeometry<GridLayoutT, std::decay_t<decltype(qty)>>::toFieldBox(box, qty,
                                                                                        layout));
            }

            return std::make_shared<FieldOverlap>(destinationBoxes, offset);
        });

        return std::make_shared<TensorFieldOverlap<rank>>(std::move(overlaps));
    }



    static SAMRAI::hier::Box toFieldBox(SAMRAI::hier::Box box, tensor_t qty,
                                        GridLayoutT const& layout)
    {
        SAMRAI::hier::IntVector lower = box.lower();
        SAMRAI::hier::IntVector upper = box.upper();
        using core::dirX, core::dirY, core::dirZ;

        auto const centerings
            = ConstArray<core::QtyCentering, dimension>(core::QtyCentering::primal);

        core::for_N<dimension>([&](auto i) {
            box.setLower(i, lower[i]);
            auto const is_primal = (centerings[i] == core::QtyCentering::primal) ? 1 : 0;
            box.setUpper(i, upper[i] + is_primal);
        });


        return box;
    }

    static GridLayoutT layoutFromBox(SAMRAI::hier::Box const& box, GridLayoutT const& layout)
    {
        std::array<std::uint32_t, dimension> nbCell;
        for (std::size_t iDim = 0; iDim < dimension; ++iDim)
        {
            nbCell[iDim] = static_cast<std::uint32_t>(box.numberCells(iDim));
        }

        return GridLayoutT(layout.meshSize(), nbCell, layout.origin());
    }


private:
    GridLayoutT layout_;
    tensor_t quantity_;


    void computeDestinationBoxes_(SAMRAI::hier::BoxContainer& destinationBoxes,
                                  TensorFieldGeometry const& sourceGeometry,
                                  SAMRAI::hier::Box const& sourceMask,
                                  SAMRAI::hier::Box const& fillBox, bool const overwriteInterior,
                                  SAMRAI::hier::Transformation const& sourceOffset,
                                  SAMRAI::hier::BoxContainer const& destinationRestrictBoxes
                                  = SAMRAI::hier::BoxContainer()) const
    {
        SAMRAI::hier::Box sourceShift = sourceGeometry.ghostTensorFieldBox_ * sourceMask;
        sourceOffset.transform(sourceShift);


        bool withGhosts = true;

        core::GridLayout sourceShiftLayout = layoutFromBox(sourceShift, sourceGeometry.layout_);
        core::GridLayout fillBoxLayout     = layoutFromBox(fillBox, layout_);

        auto const& destinationBox = this->ghostTensorFieldBox_;

        SAMRAI::hier::Box const sourceBox{toFieldBox(sourceShift, quantity_, sourceShiftLayout)};

        SAMRAI::hier::Box const fillTensorField{toFieldBox(fillBox, quantity_, fillBoxLayout)};


        SAMRAI::hier::Box const together(destinationBox * sourceBox * fillTensorField);

        if (!together.empty())
        {
            if (overwriteInterior)
            {
                destinationBoxes.push_back(together);
            }
            else
            {
                destinationBoxes.removeIntersections(together, this->interiorTensorFieldBox_);
            }
        }

        if (!destinationRestrictBoxes.empty() && !destinationBoxes.empty())
        {
            SAMRAI::hier::BoxContainer restrictBoxes;
            for (auto box = destinationRestrictBoxes.begin(); box != destinationRestrictBoxes.end();
                 ++box)
            {
                restrictBoxes.push_back(toFieldBox(*box, quantity_, layoutFromBox(*box, layout_)));
            }

            destinationBoxes.intersectBoxes(restrictBoxes);
        }
    }

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    doOverlap_(TensorFieldGeometry const& destinationGeometry,
               TensorFieldGeometry const& sourceGeometry, SAMRAI::hier::Box const& sourceMask,
               SAMRAI::hier::Box const& fillBox, bool const overwriteInterior,
               SAMRAI::hier::Transformation const& sourceOffset,
               SAMRAI::hier::BoxContainer const& destinationRestrictBoxes) const
    {
        SAMRAI::hier::BoxContainer destinationBoxes;

        destinationGeometry.computeDestinationBoxes_(destinationBoxes, sourceGeometry, sourceMask,
                                                     fillBox, overwriteInterior, sourceOffset,
                                                     destinationRestrictBoxes);

        return std::make_shared<FieldOverlap>(destinationBoxes, sourceOffset);
    }
};

} // namespace PHARE::amr


#endif
