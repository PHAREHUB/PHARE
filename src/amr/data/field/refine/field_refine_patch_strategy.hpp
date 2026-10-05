#ifndef PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP

#include "amr/data/field/field_data_traits.hpp"
#include "amr/data/tensorfield/tensor_field_data_traits.hpp"

#include "core/boundary/boundary_defs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <vector>

namespace PHARE::amr
{

/**
 * @brief Strategy for filling physical boundary conditions and customizing patch refinement.
 *
 * This class implements the SAMRAI::xfer::RefinePatchStrategy interface to
 * specify how physical boundary conditions must be enforced for patches that touch
 * the domain boundaries. Refinement customization via preprocessRefine and postProcessRefine is
 * deferred to child classes.
 *
 * Each refiner is expected to hold an instance of this class, but all these instances will point to
 * the same common boundary manager.
 *
 * @tparam ScalarOrTensorFieldDataT The data type for fields or tensor fields.
 * @tparam BoundaryManagerT Manager responsible for providing boundary condition objects.
 */
template<typename ScalarOrTensorFieldDataT, typename BoundaryManagerT>
    requires(IsFieldData<ScalarOrTensorFieldDataT> || IsTensorFieldData<ScalarOrTensorFieldDataT>)
class FieldRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    static constexpr bool is_scalar        = IsFieldData<ScalarOrTensorFieldDataT>;
    static constexpr bool is_tensor        = !is_scalar;
    static constexpr std::size_t dimension = ScalarOrTensorFieldDataT::dimension;

    using field_geometry_type    = FieldGeometrySelector<ScalarOrTensorFieldDataT, is_scalar>::type;
    using gridlayout_type        = ScalarOrTensorFieldDataT::gridlayout_type;
    using grid_type              = ScalarOrTensorFieldDataT::grid_type;
    using field_type             = grid_type::field_type;
    using physical_quantity_type = BoundaryManagerT::physical_quantity_type;
    using vectorfield_type       = core::VecField<field_type, physical_quantity_type>;
    using scalar_or_tensor_field_type
        = ScalarOrTensorFieldSelector<ScalarOrTensorFieldDataT, is_scalar>::type;
    using scalar_quantity_type = physical_quantity_type::Scalar;
    using vector_quantity_type = physical_quantity_type::Vector;

    using patch_geometry_type = SAMRAI::hier::PatchGeometry;

    using boundary_type = BoundaryManagerT::boundary_type;
    using boundary_condition_type
        = core::IFieldBoundaryCondition<scalar_or_tensor_field_type, gridlayout_type>;

    /**
     * @brief Constructor.
     * @param boundary_manager Manager handling boundary conditions.
     */
    FieldRefinePatchStrategy(BoundaryManagerT& boundaryManager)
        : boundaryManager_{boundaryManager}
        , data_id_{-1}
    {
    }

    /**
     * @brief Check that the patch data identifier is registered.
     */
    void assertIDsSet() const
    {
        assert(data_id_ >= 0 && "FieldRefinePatchStrategy: IDs must be registered before use");
    }

    /**
     * @brief Register the SAMRAI patch data identifier.
     * @param field_id Integer ID from the SAMRAI variable database.
     */
    void registerIDs(int const field_id) { data_id_ = field_id; }

    void setFillPhysicalBoundaries(bool const fill) { fillPhysicalBoundaries_ = fill; }

    /**
     * @brief Apply physical boundary conditions via SAMRAI callback.
     *
     * Iterate over patch boundaries that touch a physical domain boundary and apply the appropriate
     * PHARE boundary condition to ghost regions.
     *
     * @param patch The fine patch being refined.
     * @param fill_time Simulation time for BC application.
     * @param ghost_width_to_fill Width of ghost cell layer to be filled.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch, double const fill_time,
                                       SAMRAI::hier::IntVector const& ghost_width_to_fill) override
    {
        if (!fillPhysicalBoundaries_)
            return;

        gridlayout_type const& gridLayout = ScalarOrTensorFieldDataT::getLayout(patch, data_id_);

        assert(ghost_width_to_fill <= SAMRAI::hier::IntVector(
                   static_cast<SAMRAI::tbox::Dimension>(static_cast<int>(dimension)),
                   static_cast<int>(gridLayout.options.field_ghost_width)));

        std::shared_ptr<patch_geometry_type> patchGeom = patch.getPatchGeometry();
        assert(patchGeom && "patch has no geometry.");

        auto scalarOrTensorField = [&]() {
            if constexpr (is_scalar)
            {
                return *(&(ScalarOrTensorFieldDataT::getField(patch, data_id_)));
            }
            else
            {
                return ScalarOrTensorFieldDataT::getTensorField(patch, data_id_);
            };
        }();

        // must be retrieved to pass as argument to patchGeom->getBoundaryFillBox later
        SAMRAI::hier::Box const& patch_box = patch.getBox();

        // iterations on potential boundary codimensions in [[1, dim]]
        core::for_N<dimension>([&](auto tag) {
            constexpr auto codim = tag.value + 1;

            // find all boundaries with the current codimension
            std::vector<SAMRAI::hier::BoundaryBox> const& boundaries
                = patchGeom->getCodimensionBoundaries(static_cast<int>(codim));

            // iterate on all found boundaries of given codimension
            for (SAMRAI::hier::BoundaryBox const& bBox : boundaries)
            {
                // retrieve the localBox of ghost that must be filled
                SAMRAI::hier::Box samraiBoxToFill
                    = patchGeom->getBoundaryFillBox(bBox, patch_box, ghost_width_to_fill);
                auto localBox = gridLayout.AMRToLocal(phare_box_from<dimension>(samraiBoxToFill));

                // get location of the currently treated boundary
                auto const currentBoundaryLocation
                    = static_cast<core::CodimNBoundaryLocation<codim>>(bBox.getLocationIndex());

                // get the "master" 1-codimensional boundary that applies at the currently
                // treated boundary: for instance corner in 2D belongs to two different
                // 1-codimensional boundaries (edges), so two boundary conditions compete there.
                // The responsibility of choosing which boundary condition prevails there is on
                // the boundaryManager. If the current boundary is itself 1-codimensional, then
                // masterBoundaryLocation = currentBoundaryLocation.
                core::BoundaryLocation const masterBoundaryLocation
                    = boundaryManager_.getMasterBoundaryLocation(currentBoundaryLocation);
                auto* const masterBoundary = boundaryManager_.getBoundary(masterBoundaryLocation);
                if (!masterBoundary)
                    throw std::runtime_error("Boundary not found.");

                // get the boundary condition for the current physical quantity.
                std::shared_ptr<boundary_condition_type> bc
                    = masterBoundary->getFieldCondition(scalarOrTensorField.physicalQuantity());
                if (!bc)
                    throw std::runtime_error("Field boundary condition not found.");

                // apply the retained boundary condition, as if the current boundary was belonging
                // to the 1-codimensional master boundary; this essentially defines which Cartesian
                // direction is considered to be the normal one.
                bc->apply(scalarOrTensorField, masterBoundaryLocation, localBox, gridLayout,
                          fill_time);
            }
        });
    }


    SAMRAI::hier::IntVector
    getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector{dim, 1};
    }


    void preprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                          SAMRAI::hier::Box const& fine_box,
                          SAMRAI::hier::IntVector const& ratio) override
    {
    }


    void postprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                           SAMRAI::hier::Box const& fine_box,
                           SAMRAI::hier::IntVector const& ratio) override
    {
    }


protected:
    BoundaryManagerT& boundaryManager_; //!< a reference to the boundary manager
    int data_id_; //!< the id of the resource to which this refine patch strategy is attached
    bool fillPhysicalBoundaries_ = true;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP
