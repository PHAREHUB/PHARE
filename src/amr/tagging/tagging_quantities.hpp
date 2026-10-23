#ifndef PHARE_TAGGING_QUANTITIES_HPP
#define PHARE_TAGGING_QUANTITIES_HPP

#include "core/utilities/types.hpp"

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
 * @brief Read the names and thresholds of the tagged quantities from the config dict.
 *
 * The quantities are the entries "Q0" to "Q<nbr_quantities - 1>", each with a "name" and a
 * "threshold", as written by pharein. Without "nbr_quantities", B is tagged with the threshold
 * "threshold" (default 0.1).
 *
 * @param dict Config dict at path "simulation/AMR/refinement/tagging"
 * @return The threshold of each quantity, keyed by quantity name
 * @throws std::runtime_error if a quantity name is given more than once
 */
inline std::map<std::string, double> parseQuantities(initializer::PHAREDict const& dict)
{
    if (!dict.contains("nbr_quantities"))
        return {{"B", cppdict::get_value(dict, "threshold", 0.1)}};

    std::map<std::string, double> quantities;
    auto const nbrQuantities = dict["nbr_quantities"].template to<int>();
    for (int i = 0; i < nbrQuantities; ++i)
    {
        auto const& quantity = dict["Q" + std::to_string(i)];
        auto const name      = quantity["name"].template to<std::string>();
        if (!quantities.emplace(name, quantity["threshold"].template to<double>()).second)
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
 * @brief The largest number of components of a tagged quantity, the one of a vector field.
 *
 * @tparam Model The physical model
 */
template<typename Model>
auto constexpr max_tagged_components_v
    = std::decay_t<decltype(magneticField(std::declval<Model&>()))>::size();


/**
 * @brief A function returning the components of a quantity, taken from a model set on a patch.
 *
 * @tparam Model The physical model
 */
template<typename Model>
using quantity_accessor_t = std::vector<tagged_component_t<Model>> (*)(Model&);


/**
 * @brief A tagged quantity, resolved from its name: how to get its components from
 * the model, along with its threshold.
 *
 * @tparam Model The physical model
 */
template<typename Model>
struct TaggedQuantityAccessor
{
    quantity_accessor_t<Model> components; ///< the components of the quantity on the current patch
    double threshold;                      ///< the criterion's threshold for this quantity
};


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
    static_assert(max_tagged_components_v<Model> >= 3);

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


/**
 * @brief Resolve each quantity name to the accessor to its components, keeping its threshold.
 *
 * @tparam Model The physical model
 * @param quantities The threshold of each quantity, keyed by quantity name
 * @return the accessor and threshold of each quantity, in name order
 * @throws std::runtime_error as resolveQuantity() does
 */
template<typename Model>
std::vector<TaggedQuantityAccessor<Model>>
resolveQuantities(std::map<std::string, double> const& quantities)
{
    std::vector<TaggedQuantityAccessor<Model>> resolved;
    for (auto const& [name, threshold] : quantities)
        resolved.emplace_back(resolveQuantity<Model>(name), threshold);
    return resolved;
}

/**
 * @brief Read the tagged quantities from the config dict and resolve them.
 *
 * @tparam Model The physical model
 * @param dict config dict at path "simulation/AMR/refinement/tagging"
 * @return the accessor and threshold of each quantity, in name order
 * @throws std::runtime_error as parseQuantities() and resolveQuantity() do
 */
template<typename Model>
std::vector<TaggedQuantityAccessor<Model>> resolveQuantities(initializer::PHAREDict const& dict)
{
    return resolveQuantities<Model>(parseQuantities(dict));
}

} // namespace PHARE::amr::tagging

#endif
