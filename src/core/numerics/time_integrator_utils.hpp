#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_UTILS_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_UTILS_HPP


#include "core/utilities/index/index.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

namespace PHARE::core
{
template<typename Float, typename State>
struct RKPair
{
    Float weight;
    State& state;
};

template<typename GridLayout>
class RKUtils
{
    constexpr static auto dimension = GridLayout::dimension;


public:
    RKUtils(GridLayout const& layout)
        : layout_{layout}
    {
    }


    template<typename ReturnState>
    void operator()(ReturnState& res, auto const... pairs) const
    {
        auto result_fields = getFieldTuples_(res);

        auto weight_tuple = std::make_tuple(pairs.weight...);

        auto state_field_tuples = std::make_tuple(getFieldTuples_(pairs.state)...);

        constexpr auto num_fields = std::tuple_size_v<std::decay_t<decltype(result_fields)>>;

        for_N<num_fields>([&](auto i) {
            layout_.evalOnGhostBox(std::get<i>(result_fields), [&](auto... indices) {
                RKstep_(result_fields, weight_tuple, state_field_tuples, i, {indices...});
            });
        });
    }

private:
    template<typename State>
    static auto getFieldTuples_(State& state)
    {
        return std::forward_as_tuple(state.rho, state.rhoV(Component::X), state.rhoV(Component::Y),
                                     state.rhoV(Component::Z), state.B1(Component::X),
                                     state.B1(Component::Y), state.B1(Component::Z), state.Etot1);
    }

    template<typename ReturnState, typename WeightsTuple, typename StatesTuple, typename IndexType>
    static void RKstep_(ReturnState& res, WeightsTuple const& weights, StatesTuple const& states,
                        IndexType field_index, MeshIndex<dimension> const index)
    {
        auto sum = 0.0;

        constexpr auto num_terms = std::tuple_size_v<std::decay_t<decltype(weights)>>;

        for_N<num_terms>([&](auto j) {
            sum += std::get<j>(weights) * std::get<field_index>(std::get<j>(states))(index);
        });

        std::get<field_index>(res)(index) = sum;
    }


    GridLayout layout_;
};

} // namespace PHARE::core

#endif
