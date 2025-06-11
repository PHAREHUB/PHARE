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


#include <cassert>
#include <memory>
#include <type_traits>

namespace PHARE::amr
{


template<std::size_t dimension, std::size_t rank = 1>
class TensorFieldGeometryBase : public SAMRAI::hier::BoxGeometry
{
    using FieldGeometryBase_t = FieldGeometryBase<dimension>;

    static constexpr std::size_t N = core::detail::tensor_field_dim_from_rank<rank>();

public:
    virtual ~TensorFieldGeometryBase() {}
    TensorFieldGeometryBase(std::array<std::shared_ptr<FieldGeometryBase_t>, N>&& geoms)
        // maybe add a check that all geoms have the same patchBox?
        : patchBox{geoms[0]->patchBox}
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            components_[i] = std::move(geoms[i]);
        }
    }

    std::array<SAMRAI::hier::Box, N> interiorTensorFieldBox() const
    {
        return core::for_N<N, core::for_N_R_mode::make_array>(
            [&](auto i) { return components_[i]->interiorFieldBox(); });
    }

    SAMRAI::hier::Box const patchBox;

private:
    std::array<std::shared_ptr<FieldGeometryBase_t>, N> components_;
};

template<std::size_t rank, typename GridLayoutT, typename PhysicalQuantity>
class TensorFieldGeometry : public TensorFieldGeometryBase<GridLayoutT::dimension, rank>
{
    using tensor_t        = typename PhysicalQuantity::template TensorType<rank>;
    using FieldGeometry_t = FieldGeometry<GridLayoutT, typename PhysicalQuantity::Scalar>;

    auto static make_geoms(SAMRAI::hier::Box const& box, GridLayoutT const& layout,
                           tensor_t const qty)
    {
        auto qts         = PhysicalQuantity::componentsQuantities(qty);
        auto components_ = core::for_N<N, core::for_N_R_mode::make_array>([&](auto i) {
            return std::make_shared<FieldGeometry<GridLayoutT, std::decay_t<decltype(qts[i])>>>(
                box, layout, qts[i]);
        });

        auto base_ptr = core::for_N<N, core::for_N_R_mode::make_array>([&](auto i) {
            return std::static_pointer_cast<FieldGeometryBase<GridLayoutT::dimension>>(
                components_[i]);
        });

        return std::make_pair(std::move(base_ptr), std::move(components_));
    }

public:
    using Super                            = TensorFieldGeometryBase<GridLayoutT::dimension, rank>;
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    static constexpr std::size_t interp_order = GridLayoutT::interp_order;

    static constexpr auto N = core::detail::tensor_field_dim_from_rank<rank>();

    TensorFieldGeometry(SAMRAI::hier::Box const& box, GridLayoutT const& layout, tensor_t const qty)
        : TensorFieldGeometry(box, layout, qty, make_geoms(box, layout, qty))
    {
    }


    NO_DISCARD auto& operator[](std::size_t i) { return components_[i]; }
    NO_DISCARD auto& operator[](std::size_t i) const { return components_[i]; }


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

        auto overlaps = core::for_N<N, core::for_N_R_mode::make_array>([&](auto i) {
            auto overlap = components_[i]->calculateOverlap(
                *destinationCast[i], *sourceCast[i], sourceMask, fillBox, overwriteInterior,
                sourceOffset, retry, destinationRestrictBoxes);

            return std::dynamic_pointer_cast<FieldOverlap>(overlap);
        });

        return std::make_shared<TensorFieldOverlap<rank>>(std::move(overlaps));
    }




    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    setUpOverlap(SAMRAI::hier::BoxContainer const& boxes,
                 SAMRAI::hier::Transformation const& offset) const final
    {
        auto overlaps = core::for_N<N, core::for_N_R_mode::make_array>([&](auto i) {
            auto overlap = components_[i]->setUpOverlap(boxes, offset);
            return std::dynamic_pointer_cast<FieldOverlap>(overlap);
        });

        return std::make_shared<TensorFieldOverlap<rank>>(std::move(overlaps));
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
    // helper constructor to make sure instantiation happens in the right order
    TensorFieldGeometry(SAMRAI::hier::Box const& box, GridLayoutT const& layout, tensor_t const qty,
                        auto geoms)
        : Super(std::move(geoms.first))
        , components_(std::move(geoms.second))
    {
        for (auto component : components_)
        {
            if (!component)
            {
                throw std::runtime_error("TensorFieldGeometry: component is null");
            }
        }
    }

    std::array<std::shared_ptr<FieldGeometry_t>, N> components_;
};

} // namespace PHARE::amr


#endif
