#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/differential_operators.hpp"

#include <cstdint>

namespace PHARE::core
{


enum class AmpereMode { ShrinkedGhost, GrownPhysical };

struct AmpereBox // structural type — usable as NTTP
{
    AmpereMode mode      = AmpereMode::ShrinkedGhost;
    std::uint32_t amount = 1;
};


template<typename GridLayout, AmpereBox box = AmpereBox{}>
class Ampere
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    Ampere(GridLayout const& layout)
        : layout_{layout}
    {
    }


    template<typename VecField>
    void operator()(VecField const& B, VecField& J)
    {
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        Point<std::uint32_t, dimension> amount;

        for (size_t i = 0; i < dimension; ++i)
        {
            amount[i] = box.amount;
        }

        // evalOn* calls fn(point, trailingArgs...); trailing args (Jc, curlB component)
        // are passed through it rather than captured, so the lambdas stay capture-free
        auto const eval = [&](auto& Jc, auto&& fn, auto const& curlBc) {
            if constexpr (box.mode == AmpereMode::ShrinkedGhost)
                layout_.evalOnShrinkedGhostBox(Jc, amount, fn, Jc, curlBc);
            else
                layout_.evalOnBiggerBox(Jc, amount, fn, Jc, curlBc);
        };

        auto const curlB = curl(layout_, B);

        auto constexpr JEq
            = [](auto const& ijk, auto& Jc, auto const& curlBc) { Jc(ijk) = curlBc(ijk); };

        eval(Jx, JEq, curlB.x);
        eval(Jy, JEq, curlB.y);
        eval(Jz, JEq, curlB.z);
    }


private:
    GridLayout layout_;
};

} // namespace PHARE::core
#endif
