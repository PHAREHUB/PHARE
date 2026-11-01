#ifndef PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_UPWIND_EDGE_B_HPP
#define PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_UPWIND_EDGE_B_HPP

#include "core/def.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include <cstddef>

namespace PHARE::core
{
struct UpwindEdgeStates
{
    double aL, aR, dL, dR, BL, BR;

    NO_DISCARD double edge() const { return aL * BL + aR * BR; }
};


template<typename GridLayout, template<typename> typename Reconstruction>
struct UpwindEdgeB
{
    static constexpr auto dimension = GridLayout::dimension;

    template<Component edge, Component component>
    static UpwindEdgeStates states(auto const& ct_state, auto const& B,
                                   MeshIndex<dimension> const& idx)
    {
        constexpr auto e = static_cast<std::size_t>(edge);
        constexpr auto c = static_cast<std::size_t>(component);
        static_assert(e < 3 and c < 3 and e != c);
        constexpr auto d = 3 - e - c;
        static_assert(d < dimension);

        constexpr auto across = static_cast<Direction>(d);

        auto const [BL, BR]
            = Reconstruction<GridLayout>::template reconstruct<across>(B(component), idx);

        auto const weight = [&](auto const& w) {
            if constexpr (c < dimension)
                return 0.5
                       * (w(idx)
                          + w(GridLayout::template previous<static_cast<Direction>(c)>(idx)));
            else
                return w(idx);
        };

        return {weight(aL_<across>(ct_state)),
                weight(aR_<across>(ct_state)),
                weight(dL_<across>(ct_state)),
                weight(dR_<across>(ct_state)),
                BL,
                BR};
    }

    template<Component edge, Component component>
    static double at(auto const& ct_state, auto const& B, MeshIndex<dimension> const& idx)
    {
        return states<edge, component>(ct_state, B, idx).edge();
    }

private:
    template<Direction direction>
    static auto const& aL_(auto const& ct_state)
    {
        if constexpr (direction == Direction::X)
            return ct_state.aL_x;
        else if constexpr (direction == Direction::Y)
            return ct_state.aL_y;
        else
            return ct_state.aL_z;
    }

    template<Direction direction>
    static auto const& aR_(auto const& ct_state)
    {
        if constexpr (direction == Direction::X)
            return ct_state.aR_x;
        else if constexpr (direction == Direction::Y)
            return ct_state.aR_y;
        else
            return ct_state.aR_z;
    }

    template<Direction direction>
    static auto const& dL_(auto const& ct_state)
    {
        if constexpr (direction == Direction::X)
            return ct_state.dL_x;
        else if constexpr (direction == Direction::Y)
            return ct_state.dL_y;
        else
            return ct_state.dL_z;
    }

    template<Direction direction>
    static auto const& dR_(auto const& ct_state)
    {
        if constexpr (direction == Direction::X)
            return ct_state.dR_x;
        else if constexpr (direction == Direction::Y)
            return ct_state.dR_y;
        else
            return ct_state.dR_z;
    }
};

} // namespace PHARE::core

#endif
