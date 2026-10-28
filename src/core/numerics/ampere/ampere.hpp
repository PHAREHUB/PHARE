#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"
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
    void operator()(VecField const& B, VecField& J) _PHARE_ALL_FN_
    {
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        Point<std::uint32_t, dimension> const amount{
            ConstArray<std::uint32_t, dimension>(box.amount)};

        auto const eval = [&] _PHARE_ALL_FN_(auto& Jc, auto&& fn) {
            if constexpr (box.mode == AmpereMode::ShrinkedGhost)
                layout_.evalOnShrinkedGhostBox(Jc, amount, fn);
            else
                layout_.evalOnBiggerBox(Jc, amount, fn);
        };

        using Field_t    = std::remove_cvref_t<decltype(Jx)>::Super;
        auto const curlB = curl(layout_, B);
        Field_t jx = Jx, jy = Jy, jz = Jz;

        eval(Jx,
             [=] _PHARE_ALL_FN_(auto const&... ijk) mutable { jx(ijk...) = curlB.x({ijk...}); });
        eval(Jy,
             [=] _PHARE_ALL_FN_(auto const&... ijk) mutable { jy(ijk...) = curlB.y({ijk...}); });
        eval(Jz,
             [=] _PHARE_ALL_FN_(auto const&... ijk) mutable { jz(ijk...) = curlB.z({ijk...}); });
    }


private:
    GridLayout layout_;
};

} // namespace PHARE::core

#endif
