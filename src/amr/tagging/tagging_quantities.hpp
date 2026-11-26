#ifndef PHARE_TAGGING_QUANTITIES_HPP
#define PHARE_TAGGING_QUANTITIES_HPP

#include "core/utilities/meta/type_meta.hpp"

#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"

#include "initializer/data_provider.hpp"

#include <map>
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
/**
 * @brief The threshold of a tagged quantity, along with its config dict, from which the criterion
 * reads its per-quantity parameters.
 */
struct QuantitySpec
{
    double threshold;                  ///< the criterion's threshold for this quantity
    initializer::PHAREDict parameters; ///< the quantity's config dict, empty for the fallback on B
};

/**
 * @brief Read the names, thresholds and config dicts of the tagged quantities from the config dict.
 *
 * The quantities are the entries "Q0" to "Q<nbr_quantities - 1>", each with a "name", a
 * "threshold" and the per-quantity parameters of the criterion, as written by pharein. Without
 * "nbr_quantities", B is tagged with the threshold "threshold" (default 0.1) and no parameter.
 *
 * @param dict Config dict at path "simulation/AMR/refinement/tagging"
 * @return The threshold and config dict of each quantity, keyed by quantity name
 * @throws std::runtime_error if a quantity name is given more than once
 */
inline std::map<std::string, QuantitySpec> parseQuantities(initializer::PHAREDict const& dict)
{
    if (!dict.contains("nbr_quantities"))
        return {{"B", QuantitySpec{cppdict::get_value(dict, "threshold", 0.1), {}}}};

    std::map<std::string, QuantitySpec> quantities;
    auto const nbrQuantities = dict["nbr_quantities"].template to<int>();
    for (int i = 0; i < nbrQuantities; ++i)
    {
        auto const& quantity = dict["Q" + std::to_string(i)];
        auto const name      = quantity["name"].template to<std::string>();
        QuantitySpec spec{quantity["threshold"].template to<double>(), quantity};
        if (!quantities.emplace(name, std::move(spec)).second)
            throw std::runtime_error("tagging quantity '" + name + "' is given more than once");
    }
    return quantities;
}


/**
 * @brief The magnetic field of the model.
 *
 * @tparam Model The physical model
 * @param model The physical model
 * @return a reference to the model's magnetic field
 */
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


/**
 * @brief The field type of the tagged components, the one of a magnetic field component.
 *
 * @tparam Model The physical model
 */
template<typename Model>
using tagged_field_t = std::decay_t<
    decltype(magneticField(std::declval<Model&>()).getComponent(core::Component::X))>;

/**
 * @brief A tagged component: a copy of the component's field view (data pointer and shape), only
 * valid on the patch the model was set on when it was taken.
 *
 * @tparam Model The physical model
 */
template<typename Model>
using tagged_component_t = tagged_field_t<Model>;


/**
 * @brief A function returning the components of a quantity, taken from a model set on a patch.
 *
 * @tparam Model The physical model
 */
template<typename Model>
using quantity_accessor_t = std::vector<tagged_component_t<Model>> (*)(Model&);


/**
 * @brief The quantities available for tagging in the model, each name mapped to the accessor of its
 * components.
 *
 * B, Bx, By and Bz for both models, plus rho (or density), rhoV (or momentum), rhoVx, rhoVy, rhoVz
 * and Etot for the MHD model.
 *
 * @tparam Model The physical model
 * @return the list of (name, accessor) pairs
 */
template<typename Model>
auto taggableQuantities()
{
    using Tagged     = tagged_component_t<Model>;
    using Components = std::vector<Tagged>;
    using Accessor   = quantity_accessor_t<Model>;
    using core::Component;

    static constexpr auto scalar = [](auto const& field) { return Components{Tagged{field}}; };
    static constexpr auto one    = [](auto const& vecfield, Component c) {
        return Components{Tagged{vecfield.getComponent(c)}};
    };
    static constexpr auto all = [](auto const& vecfield) {
        return Components{Tagged{vecfield.getComponent(Component::X)},
                          Tagged{vecfield.getComponent(Component::Y)},
                          Tagged{vecfield.getComponent(Component::Z)}};
    };

    std::vector<std::pair<std::string, Accessor>> quantities{
        {"B", [](Model& model) { return all(magneticField(model)); }},
        {"Bx", [](Model& model) { return one(magneticField(model), Component::X); }},
        {"By", [](Model& model) { return one(magneticField(model), Component::Y); }},
        {"Bz", [](Model& model) { return one(magneticField(model), Component::Z); }}};

    if constexpr (solver::is_mhd_model_v<Model>)
    {
        quantities.insert(
            quantities.end(),
            {{"rho", [](Model& model) { return scalar(model.state.rho); }},
             {"density", [](Model& model) { return scalar(model.state.rho); }},
             {"rhoV", [](Model& model) { return all(model.state.rhoV); }},
             {"momentum", [](Model& model) { return all(model.state.rhoV); }},
             {"rhoVx", [](Model& model) { return one(model.state.rhoV, Component::X); }},
             {"rhoVy", [](Model& model) { return one(model.state.rhoV, Component::Y); }},
             {"rhoVz", [](Model& model) { return one(model.state.rhoV, Component::Z); }},
             {"Etot", [](Model& model) { return scalar(model.state.Etot); }}});
    }

    return quantities;
}


/**
 * @brief The names of quantities of the model that are not available for tagging yet, so that they
 * get a clearer error than unknown names.
 *
 * @tparam Model The physical model
 * @return the list of names
 */
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


/**
 * @brief The accessor to the components of the quantity of the given name.
 *
 * @tparam Model The physical model
 * @param name The name of the quantity
 * @return the accessor to the quantity's components
 * @throws std::runtime_error if the quantity is not available for tagging yet, or unknown for this
 * model (the message then lists the allowed names)
 */
template<typename Model>
quantity_accessor_t<Model> resolveQuantity(std::string const& name)
{
    auto const quantities = taggableQuantities<Model>();

    auto const found = std::find_if(std::begin(quantities), std::end(quantities),
                                    [&](auto const& quantity) { return quantity.first == name; });
    if (found != std::end(quantities))
        return found->second;

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


} // namespace PHARE::amr::tagging

#endif
