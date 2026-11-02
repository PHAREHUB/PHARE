#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_EXTERNAL_FIELD_ON_LAYOUT_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_EXTERNAL_FIELD_ON_LAYOUT_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/models/external_field_updater.hpp"
#include "core/utilities/types.hpp"

#include <array>
#include <concepts>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>

namespace PHARE::core
{
/**
 * @brief Evaluates the external field B0 on any grid layout, in memory it owns.
 *
 * Boundary conditions run on patches that may carry no external field patch data (the
 * temporary patches SAMRAI creates while regridding), so they recompute B0 from the model's
 * updater on the patch layout at the fill time. This gives the same values as the B0 stored by
 * the model, since it is the same discrete curl of the same potential on the same layout.
 *
 * @tparam VecFieldT vecfield implementation
 * @tparam GridLayoutT grid layout implementation
 */
template<typename VecFieldT, typename GridLayoutT>
class ExternalFieldOnLayout
{
public:
    using updater_type = IExternalFieldUpdater<VecFieldT, GridLayoutT>;
    using field_type   = VecFieldT::field_type;
    using value_type   = field_type::value_type;
    using tensor_type  = VecFieldT::tensor_t;

    static constexpr std::size_t dimension = GridLayoutT::dimension;

    using grid_type
        = Grid<NdArrayVector<dimension, value_type>, typename field_type::physical_quantity_type>;
    static_assert(std::same_as<typename grid_type::field_type, field_type>);

    explicit ExternalFieldOnLayout(std::shared_ptr<updater_type> updater)
        : updater_{std::move(updater)}
    {
        if (!updater_)
            throw std::runtime_error("ExternalFieldOnLayout: null external field updater");
    }

    /**
     * @brief compute B0 on the full ghost box of a layout, at a given time
     *
     * @return a view on the computed B0, valid until the next call
     */
    VecFieldT const& B0(GridLayoutT const& layout, double const time)
    {
        allocate_(B0_, B0grids_, layout);
        updater_->computeB0(B0_, layout, time);
        return B0_;
    }

private:
    using grids_type = std::array<grid_type, VecFieldT::size()>;

    static void allocate_(VecFieldT& vecfield, std::optional<grids_type>& grids,
                          GridLayoutT const& layout)
    {
        grids.emplace(for_N<VecFieldT::size(), for_N_R_mode::make_array>([&](auto i) {
            auto const qty = vecfield[i].physicalQuantity();
            return grid_type{vecfield[i].name(), qty, layout.allocSize(qty)};
        }));
        vecfield.setBuffer(&*grids);
    }

    std::shared_ptr<updater_type> updater_;
    VecFieldT B0_{"boundary_external_field_B0", tensor_type::B};
    std::optional<grids_type> B0grids_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_EXTERNAL_FIELD_ON_LAYOUT_HPP
