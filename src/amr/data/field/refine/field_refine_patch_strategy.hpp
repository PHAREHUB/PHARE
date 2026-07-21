#ifndef PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP

#include "amr/data/field/field_data.hpp"
#include "amr/data/field/field_data_traits.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "amr/data/tensorfield/tensor_field_data_traits.hpp"

#include "core/boundary/boundary_defs.hpp"
#include "core/data/patch_field_accessor.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/boundary_condition_context.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <unordered_map>

namespace PHARE::amr
{

/**
 * @brief Strategy for filling physical boundary conditions and customizing patch refinment.
 *
 * This class implements the SAMRAI::xfer::RefinePatchStrategy interface to
 * specify how physical boundary conditions must be enforced for patches that touch
 * the domain boundaries. Refinement customization is deferred to child classes.
 *
 * @tparam ResMan The resources manager type.
 * @tparam ScalarOrTensorFieldDataT The data type for fields or tensor fields.
 * @tparam BoundaryManagerT Manager responsible for providing boundary condition objects.
 */
template<typename ResMan, typename ScalarOrTensorFieldDataT, typename BoundaryManagerT>
    requires(IsFieldData<ScalarOrTensorFieldDataT> || IsTensorFieldData<ScalarOrTensorFieldDataT>)
class FieldRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    static constexpr bool is_scalar   = IsFieldData<ScalarOrTensorFieldDataT>;
    static constexpr bool is_tensor   = !is_scalar;
    static constexpr size_t dimension = ScalarOrTensorFieldDataT::dimension;

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

    using patch_geometry_type           = SAMRAI::hier::PatchGeometry;
    using cartesian_patch_geometry_type = SAMRAI::geom::CartesianPatchGeometry;

    using boundary_type = BoundaryManagerT::boundary_type;
    using boundary_condition_type
        = core::IFieldBoundaryCondition<scalar_or_tensor_field_type, gridlayout_type>;
    using scalar_id_map_type     = std::unordered_map<scalar_quantity_type, int>;
    using vector_id_map_type     = std::unordered_map<vector_quantity_type, int>;
    using scalar_field_data_type = FieldData<gridlayout_type, grid_type, scalar_quantity_type>;
    using vector_field_data_type
        = TensorFieldData<1, gridlayout_type, grid_type, physical_quantity_type>;

    /**
     * @brief Concrete accessor to retrieve any field from a SAMRAI patch by physical quantity.
     *
     * Implements the core::IPatchFieldAccessor interface. Constructed once per
     * setPhysicalBoundaryConditions call and passed to boundary condition apply() methods,
     * allowing coupled BCs (e.g. inflow/outflow) to read other fields.
     *
     * Defined as a nested class to avoid heavy external template parameters.
     */
    class PatchFieldAccessor : public core::IPatchFieldAccessor<field_type, physical_quantity_type>
    {
    public:
        PatchFieldAccessor(SAMRAI::hier::Patch const& patch, scalar_id_map_type const& scalarIds,
                           vector_id_map_type const& vectorIds)
            : patch_{patch}
            , scalarIds_{scalarIds}
            , vectorIds_{vectorIds}
        {
        }

        field_type& getField(scalar_quantity_type qty) const override
        {
            auto it = scalarIds_.find(qty);
            if (it == scalarIds_.end())
                throw core::PatchFieldAccessorError(
                    "PatchFieldAccessor: scalar quantity not registered");
            // SAMRAI runs setPhysicalBoundaryConditions on temporary single-quantity
            // interpolation patches too: a sibling read by a coupled BC is registered
            // (in the id-map) but not allocated there. Surface that as the accessor error
            // the strategy catches to fall back to a sibling-free Neumann fill.
            if (!patch_.checkAllocated(it->second))
                throw core::PatchFieldAccessorError(
                    "PatchFieldAccessor: scalar quantity not allocated on patch");
            return *(&(scalar_field_data_type::getField(patch_, it->second)));
        }

        vectorfield_type getVecField(vector_quantity_type qty) const override
        {
            auto it = vectorIds_.find(qty);
            if (it == vectorIds_.end())
                throw core::PatchFieldAccessorError(
                    "PatchFieldAccessor: vector quantity not registered");
            if (!patch_.checkAllocated(it->second))
                throw core::PatchFieldAccessorError(
                    "PatchFieldAccessor: vector quantity not allocated on patch");
            return vector_field_data_type::getTensorField(patch_, it->second);
        }

    private:
        SAMRAI::hier::Patch const& patch_;
        scalar_id_map_type const& scalarIds_;
        vector_id_map_type const& vectorIds_;
    };

    using patch_field_accessor_type = PatchFieldAccessor;

    /**
     * @brief Constructor.
     * @param resources_manager Simulation resources manager.
     * @param boundary_manager Manager handling boundary conditions.
     */
    FieldRefinePatchStrategy(ResMan& resourcesManager, BoundaryManagerT& boundaryManager)
        : rm_{resourcesManager}
        , boundaryManager_{boundaryManager}
        , data_id_{-1}
        , all_scalar_ids_{}
        , all_vector_ids_{}
        , old_scalar_ids_{}
        , old_vector_ids_{}
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
     * @param all_scalar_ids id-map of scalar fields exposed to BC appliers as the *current* state.
     * @param all_vector_ids id-map of vector fields exposed to BC appliers as the *current* state.
     * @param old_scalar_ids id-map of scalar fields exposed as the *previous-substage* state.
     * @param old_vector_ids id-map of vector fields exposed as the *previous-substage* state.
     */
    void registerIDs(int const field_id, scalar_id_map_type all_scalar_ids = {},
                     vector_id_map_type all_vector_ids = {}, scalar_id_map_type old_scalar_ids = {},
                     vector_id_map_type old_vector_ids = {})
    {
        data_id_        = field_id;
        all_scalar_ids_ = std::move(all_scalar_ids);
        all_vector_ids_ = std::move(all_vector_ids);
        old_scalar_ids_ = std::move(old_scalar_ids);
        old_vector_ids_ = std::move(old_vector_ids);
    }

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
    void
    setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch, double const fill_time,
                                  SAMRAI::hier::IntVector const& /*ghost_width_to_fill*/) override
    {
        gridlayout_type const& gridLayout = ScalarOrTensorFieldDataT::getLayout(patch, data_id_);

        /// @todo SAMRAI does not pass a `ghost_width_to_fill` consistent with
        /// `gridLayout.nbrGhosts()` beyond L0, so we ignore the argument and refill the whole
        /// ghost layer. Making it consistent would require overriding getRefineOpStencilWidth,
        /// deferred to avoid perturbing the always-return-1 interpolation stencil.
        SAMRAI::hier::IntVector const ghost_width_to_fill{
            static_cast<SAMRAI::tbox::Dimension>(static_cast<int>(dimension)),
            static_cast<int>(gridLayout.nbrGhosts())};

        // no check this is a valid cast
        std::shared_ptr<cartesian_patch_geometry_type> patchGeom
            = std::static_pointer_cast<cartesian_patch_geometry_type>(patch.getPatchGeometry());

        // `*(&(getField(...)))` extracts the non-owning Field view out of the patch-owned
        // Grid via Grid::operator&() (returns &field_). `auto` then copies that lightweight
        // view, which still aliases the patch buffer, so bc->apply(...) writes to patch data.
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

        // build two accessors: one for the current substage state, one for the previous one.
        // State-aware BCs (e.g. NSCBC / LODI characteristic outflow) read from `accessor_old`,
        // integrate over `dt`, and write into ghost cells reachable via `accessor_new`.
        patch_field_accessor_type fieldAccessor{patch, all_scalar_ids_, all_vector_ids_};
        patch_field_accessor_type fieldAccessorOld{patch, old_scalar_ids_, old_vector_ids_};
        core::BoundaryConditionContext<field_type, physical_quantity_type> const ctx{
            fieldAccessor, fieldAccessorOld, fill_time};

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

                // get the primary 1-codimensional boundary that applies at the currently treated
                // boundary. If the current boundary is itself 1-codimensional, then
                // masterBoundaryLocation = currentBoundaryLocation
                core::BoundaryLocation const masterBoundaryLocation
                    = boundaryManager_.getMasterBoundaryLocation(currentBoundaryLocation);
                std::shared_ptr<boundary_type> masterBoundary
                    = boundaryManager_.getBoundary(masterBoundaryLocation);
                if (!masterBoundary)
                    throw std::runtime_error("Boundary not found.");

                // get the boundary condition for the current physical quantity. The same normal
                // condition is applied during the regular advance and while a fine level is
                // (re)filled at regrid/init: the value-prescribed inflow conditions (e.g. the
                // divergence-free transverse Dirichlet B) produce valid outside-domain ghosts
                // without needing the fine interior, and the extrapolating ones (e.g. open/outflow
                // B) read the freshly (re)filled interior.
                std::shared_ptr<boundary_condition_type> bc
                    = masterBoundary->getFieldCondition(scalarOrTensorField.physicalQuantity());
                if (!bc)
                    throw std::runtime_error("Field boundary condition not found.");

                // apply the boundary condition as if the current boundary was belonging to the
                // primary boundary.
                //
                // SAMRAI invokes this callback not only on the real level patches (which carry the
                // full MHD state) but also on temporary, single-quantity patches it builds for
                // cross-level (coarse->fine) interpolation. PHARE's coupled MHD conditions read
                // sibling fields off the patch (energy BC: rho/P/rhoV/B); those
                // siblings are not allocated on the interpolation temp patches and the accessor
                // throws there. The temp-patch ghosts still feed the fine level via interpolation,
                // so they must not be left at the NaN sentinel: fall back to a sibling-free
                // zero-gradient (Neumann) fill for that quantity. The real level patches, which
                // carry the full state, still receive the exact coupled condition.
                //
                // Only the missing-sibling case (PatchFieldAccessorError) is caught here; any other
                // exception is a genuine fault and must propagate rather than be masked as a silent
                // Neumann fill.
                try
                {
                    bc->apply(scalarOrTensorField, masterBoundaryLocation, localBox, gridLayout,
                              ctx);
                }
                catch (core::PatchFieldAccessorError const& e)
                {
                    PHARE_LOG_LINE_SS(
                        "Neumann fallback triggered in setPhysicalBoundaryConditions"
                        << " | field=" << scalarOrTensorField.name() << " | quantity="
                        << static_cast<int>(scalarOrTensorField.physicalQuantity())
                        << " | codim=" << static_cast<int>(codim) << " | currentBoundaryLocation="
                        << static_cast<int>(currentBoundaryLocation)
                        << " | masterBoundaryLocation=" << static_cast<int>(masterBoundaryLocation)
                        << " | fill_time=" << fill_time << " | patch_box=" << patch_box
                        << " | localBox=" << localBox << " | reason=" << e.what());
                    core::FieldNeumannBoundaryCondition<scalar_or_tensor_field_type,
                                                        gridlayout_type>
                        neumannFallback;
                    neumannFallback.apply(scalarOrTensorField, masterBoundaryLocation, localBox,
                                          gridLayout, ctx);
                }
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


    static auto isNewFineFace(auto const& amrIdx, auto const dir) {}


protected:
    ResMan& rm_;
    BoundaryManagerT& boundaryManager_;
    int data_id_;
    scalar_id_map_type all_scalar_ids_;
    vector_id_map_type all_vector_ids_;
    scalar_id_map_type old_scalar_ids_;
    vector_id_map_type old_vector_ids_;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_FIELD_REFINE_PATCH_STRATEGY_HPP
