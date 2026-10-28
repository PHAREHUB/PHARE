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
    using external_field_type = core::ExternalField<VecFieldT>;
    using updater_type        = core::IExternalFieldUpdater<VecFieldT, GridLayoutT>;
    using factory_type        = core::ExternalFieldUpdaterFactory<VecFieldT, GridLayoutT>;

    ExternalFieldManager(std::string const& name, initializer::PHAREDict const& dict)
        : field_{name, dict.contains("external_field")}
        , updater_{field_.isEnabled() ? factory_type::create(dict["external_field"]) : nullptr}
    {
    }

    NO_DISCARD bool isEnabled() const { return field_.isEnabled(); }

    NO_DISCARD external_field_type& field() { return field_; }
    NO_DISCARD external_field_type const& field() const { return field_; }

    void initialize(auto& resourcesManager, auto& level, VecFieldT& scratch, double const time)
    {
        if (!isEnabled())
            return;

        for (auto const& patch : resourcesManager.enumerate(level, field_, scratch))
        {
            auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto potential    = core::view_as(scratch, VecFieldT::tensor_t::E, layout);
            (*updater_)(field_, potential, layout, time);
        }
    }

    void update(auto& resourcesManager, auto& level, VecFieldT& scratch, double const time)
    {
        if (isEnabled() and updater_->isTimeDependent())
            initialize(resourcesManager, level, scratch, time);
    }

private:
    external_field_type field_;
    std::unique_ptr<updater_type> updater_;
};

} // namespace PHARE::solver

#endif // PHARE_AMR_PHYSICAL_MODELS_EXTERNAL_FIELD_MANAGER_HPP
