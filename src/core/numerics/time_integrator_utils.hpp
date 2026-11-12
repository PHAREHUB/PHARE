#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_UTILS_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_UTILS_HPP


#include "core/utilities/algorithm.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/field/field_tiles.hpp"
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

        // Combine on the physical box only. RKUtils is an interior combiner; the ghost box
        // is filled separately by the messenger's coarse-fine time interpolation at the
        // node's abscissa (blending the endpoint ghosts here would use scheme weights, not
        // the node's time-interpolation alpha).
        for_N<num_fields>([&](auto i) {
            auto const combine = [&](auto const& layout, auto&& r, auto const&... vals) {
                auto const vals_tuple = std::forward_as_tuple(vals...);
                layout.evalOnBox(r, [&](auto... indices) {
                    MeshIndex<dimension> const index{indices...};
                    double sum = 0.0;
                    for_N<sizeof...(vals)>([&](auto j) {
                        sum += std::get<j>(weight_tuple) * std::get<j>(vals_tuple)(index);
                    });
                    r(index) = sum;
                });
            };
            std::apply(
                [&](auto const&... states) {
                    auto& res = std::get<i>(result_fields);
                    if constexpr (is_field_tile_set_v<std::decay_t<decltype(res)>>)
                        tile_exec_with_layout(combine, res, std::get<i>(states)...);
                    else
                        combine(layout_, res, std::get<i>(states)...);
                },
                state_field_tuples);
        });
    }

private:
    template<typename State>
    static auto getFieldTuples_(State& state)
    {
        return std::forward_as_tuple(state.rho, state.rhoV(Component::X), state.rhoV(Component::Y),
                                     state.rhoV(Component::Z), state.B(Component::X),
                                     state.B(Component::Y), state.B(Component::Z), state.Etot);
    }


    GridLayout layout_;
};

} // namespace PHARE::core

#endif
