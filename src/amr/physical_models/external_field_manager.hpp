#ifndef PHARE_AMR_PHYSICAL_MODELS_EXTERNAL_FIELD_MANAGER_HPP
#define PHARE_AMR_PHYSICAL_MODELS_EXTERNAL_FIELD_MANAGER_HPP

#include "amr/resources_manager/amr_utils.hpp"

#include "core/data/tensorfield/tensorfield.hpp"
#include "core/def.hpp"
#include "core/models/external_field.hpp"
#include "core/models/external_field_updater.hpp"
#include "core/models/external_field_updater_factory.hpp"

#include "initializer/data_provider.hpp"

#include <memory>
#include <string>

namespace PHARE::solver
{
template<typename VecFieldT, typename GridLayoutT>
class ExternalFieldManager
{
public:
    using state_type   = core::ExternalFieldState<VecFieldT>;
    using updater_type = core::IExternalFieldUpdater<VecFieldT, GridLayoutT>;
    using factory_type = core::ExternalFieldUpdaterFactory<VecFieldT, GridLayoutT>;

    ExternalFieldManager(std::string const& name, initializer::PHAREDict const& dict)
        : updater_{dict.contains("external_field") ? factory_type::create(dict["external_field"])
                                                   : nullptr}
        , state_{name, updater_ != nullptr, updater_ and updater_->isTimeDependent()}
    {
    }

    NO_DISCARD bool isEnabled() const { return state_.isEnabled(); }

    NO_DISCARD bool isTimeDependent() const { return state_.isTimeDependent(); }

    NO_DISCARD state_type& state() { return state_; }
    NO_DISCARD state_type const& state() const { return state_; }

    void initialize(auto& resourcesManager, auto& level, VecFieldT& scratch, double const time)
    {
        if (!isEnabled())
            return;

        for (auto const& patch : resourcesManager.enumerate(level, state_, scratch))
        {
            auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto potential    = core::view_as(scratch, VecFieldT::tensor_t::E, layout);
            (*updater_)(state_, potential, layout, time);
        }
    }

    void update(auto& resourcesManager, auto& level, VecFieldT& scratch, double const time)
    {
        if (isTimeDependent())
            initialize(resourcesManager, level, scratch, time);
    }

private:
    std::unique_ptr<updater_type> updater_;
    state_type state_;
};

} // namespace PHARE::solver

#endif // PHARE_AMR_PHYSICAL_MODELS_EXTERNAL_FIELD_MANAGER_HPP
