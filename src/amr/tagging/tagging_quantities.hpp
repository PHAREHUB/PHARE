#ifndef PHARE_TAGGING_QUANTITIES_HPP
#define PHARE_TAGGING_QUANTITIES_HPP

#include "core/utilities/types.hpp"

#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"

#include "initializer/data_provider.hpp"

#include <array>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <string_view>
#include <type_traits>

namespace PHARE::amr::tagging
{
struct QuantityThreshold
{
    std::string name;
    double threshold;

    bool operator==(QuantityThreshold const&) const = default;
};


inline std::vector<QuantityThreshold> parseQuantities(initializer::PHAREDict const& dict)
{
    if (!dict.contains("nbr_quantities"))
        return {{"B", cppdict::get_value(dict, "threshold", 0.1)}};

    std::vector<QuantityThreshold> quantities;
    auto const nbrQuantities = dict["nbr_quantities"].template to<int>();
    for (int i = 0; i < nbrQuantities; ++i)
    {
        auto const& quantity = dict["Q" + std::to_string(i)];
        quantities.emplace_back(quantity["name"], quantity["threshold"]);
    }
    return quantities;
}


template<typename Model>
auto& magneticField(Model& model)
{
    if constexpr (solver::is_hybrid_model_v<Model>)
        return model.state.electromag.B;
    else if constexpr (solver::is_mhd_model_v<Model>)
        return model.state.B;
    else
        static_assert(core::dependent_false_v<Model>);
}


template<typename Model>
using tagged_field_t = std::decay_t<
    decltype(magneticField(std::declval<Model&>()).getComponent(core::Component::X))>;

template<typename Model>
using tagged_component_t = tagged_field_t<Model>;

template<typename Model>
auto constexpr max_tagged_components_v
    = std::decay_t<decltype(magneticField(std::declval<Model&>()))>::size();


template<typename Model>
auto taggableQuantities(Model& model)
{
    using Tagged     = tagged_component_t<Model>;
    using Components = std::vector<Tagged>;
    using core::Component;

    auto const scalar = [](auto const& field) { return Components{Tagged{field}}; };
    auto const one    = [](auto const& vecfield, Component c) {
        return Components{Tagged{vecfield.getComponent(c)}};
    };
    auto const all = [](auto const& vecfield) {
        return Components{Tagged{vecfield.getComponent(Component::X)},
                          Tagged{vecfield.getComponent(Component::Y)},
                          Tagged{vecfield.getComponent(Component::Z)}};
    };

    auto const& B = magneticField(model);
    std::vector<std::pair<std::string, Components>> quantities{{"B", all(B)},
                                                               {"Bx", one(B, Component::X)},
                                                               {"By", one(B, Component::Y)},
                                                               {"Bz", one(B, Component::Z)}};

    if constexpr (solver::is_mhd_model_v<Model>)
    {
        auto const& state = model.state;
        quantities.insert(quantities.end(), {{"rho", scalar(state.rho)},
                                             {"density", scalar(state.rho)},
                                             {"rhoV", all(state.rhoV)},
                                             {"momentum", all(state.rhoV)},
                                             {"rhoVx", one(state.rhoV, Component::X)},
                                             {"rhoVy", one(state.rhoV, Component::Y)},
                                             {"rhoVz", one(state.rhoV, Component::Z)},
                                             {"Etot", scalar(state.Etot)}});
    }

    return quantities;
}


template<typename Model>
constexpr auto notAvailableQuantities()
{
    if constexpr (solver::is_hybrid_model_v<Model>)
        return std::array<std::string_view, 15>{"J",
                                                "Jx",
                                                "Jy",
                                                "Jz",
                                                "E",
                                                "Ex",
                                                "Ey",
                                                "Ez",
                                                "density",
                                                "massDensity",
                                                "chargeDensity",
                                                "velocity",
                                                "bulkVel",
                                                "momentum",
                                                "pressure"};
    else
        return std::array<std::string_view, 15>{"J",  "Jx", "Jy",       "Jz", "E",
                                                "Ex", "Ey", "Ez",       "V",  "Vx",
                                                "Vy", "Vz", "velocity", "P",  "pressure"};
}


template<typename Model>
std::vector<tagged_component_t<Model>> resolveQuantity(Model& model, std::string const& name)
{
    auto quantities = taggableQuantities(model);

    auto const found = std::find_if(std::begin(quantities), std::end(quantities),
                                    [&](auto const& quantity) { return quantity.first == name; });
    if (found != std::end(quantities))
        return std::move(found->second);

    auto constexpr notAvailable = notAvailableQuantities<Model>();
    if (std::find(std::begin(notAvailable), std::end(notAvailable), name) != std::end(notAvailable))
        throw std::runtime_error("tagging quantity '" + name
                                 + "' is not available for tagging yet");

    std::string allowed;
    for (auto const& [allowedName, _] : quantities)
        allowed += (allowed.empty() ? "" : ", ") + allowedName;
    throw std::runtime_error("tagging quantity '" + name
                             + "' is unknown for this model, allowed: " + allowed);
}


template<typename Model>
auto resolveQuantities(Model& model, std::vector<QuantityThreshold> const& quantities)
{
    std::vector<std::pair<std::vector<tagged_component_t<Model>>, double>> resolved;
    for (auto const& [name, threshold] : quantities)
    {
        resolved.emplace_back(resolveQuantity(model, name), threshold);
        if (resolved.back().first.size() > max_tagged_components_v<Model>)
            throw std::runtime_error("tagging quantity '" + name + "' has more than "
                                     + std::to_string(max_tagged_components_v<Model>)
                                     + " components");
    }
    return resolved;
}

} // namespace PHARE::amr::tagging

#endif
