#ifndef PHARE_TEST_CORE_DATA_ELECTROMAG_ELECTROMAG_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_ELECTROMAG_ELECTROMAG_FIXTURES_HPP

// #include "phare_core.hpp"

#include "core/data/electromag/electromag.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/vecfield/vecfield.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include <cmath>
#include <cassert>


namespace PHARE::core
{


template<auto layout_mde, typename EM, typename GridLayout>
auto default_em_init(EM& em, GridLayout const& patch_layout)
{
    auto constexpr static dim = GridLayout::dimension;

    auto const setter = [](auto& v, auto& layout) {
        auto const box = layout.ghostBoxFor(v);
        // PHARE_LOG_LINE_SS(box);

        for (std::size_t i = 0; i < box.upper[0]; i += 2)
        {
            if constexpr (dim == 1)
            {
                v(i) += .5;
                v(i + 1) -= .5;
            }
            else
            {
                for (std::size_t j = 0; j < box.upper[1]; j += 2)
                {
                    if constexpr (dim == 2)
                    {
                        v(i, j) += .5;
                        v(i + 1, j) -= .2;
                        v(i, j + 1) += .3;
                        v(i + 1, j + 1) -= .2;
                    }
                    else
                    {
                        if constexpr (dim == 3)
                        {
                            for (std::size_t k = 0; k < box.upper[2]; k += 2)
                            {
                                v(i, j, k) += .5;
                                v(i, j + 1, k) -= .2;
                                v(i, j, k + 1) += .5;
                                v(i, j + 1, k + 1) -= .2;
                                v(i + 1, j + 1, k + 1) -= .5;
                            }
                        }
                    }
                }
            }
        }
    };

    using enum LayoutMode;
    if constexpr (is_tiled(layout_mde))
    {
        for (auto& xyz : em.E)
            for (auto& tile : xyz())
                setter(tile, tile.layout());
        for (auto& xyz : em.B)
            for (auto& tile : xyz())
                setter(tile, tile.layout());
    }
    else
    {
        for (auto& xyz : em.E)
            setter(xyz, patch_layout);
        for (auto& xyz : em.B)
            setter(xyz, patch_layout);
    }
}

template<auto opts>
struct UsableElectromagImpl
{
    using Options      = decltype(opts);
    using GridLayout_t = Options::GridLayout_t;
    using Field_t      = Options::Field_t;
    using Super_t      = Electromag<TensorField<Field_t, HybridQuantity, /*rank*/ 1>>;
};


template<auto opts>
class UsableElectromag : public UsableElectromagImpl<opts>::Super_t
{
    using Options     = decltype(opts);
    using GridLayout_ = Options::GridLayout_t;

    void _set()
    {
        E.set_on(Super::E);
        B.set_on(Super::B);
        assert(Super::isUsable());
    }

public:
    using Super = UsableElectromagImpl<opts>::Super_t;

    template<typename GridLayout>
    UsableElectromag(GridLayout const& layout, bool const init = true)
        : Super{"EM"}
        , E{"EM_E", layout, HybridQuantity::Vector::E, 1}
        , B{"EM_B", layout, HybridQuantity::Vector::B, 1}
    {
        _set();
        if (init)
            default_em_init<Options::opts.layout_mode>(*this, layout);
    }


    template<typename GridLayout>
    UsableElectromag(GridLayout const& layout, initializer::PHAREDict const& dict)
        : Super{dict}
        , E{"EM_E", layout, HybridQuantity::Vector::E, 0}
        , B{"EM_B", layout, HybridQuantity::Vector::B, 0}
    {
        _set();
    }


    UsableElectromag(UsableElectromag&& that)
        : Super{std::forward<Super>(that)}
        , E{std::move(that.E)}
        , B{std::move(that.B)}
    {
        _set();
    }


    UsableElectromag(Super const&) = delete;
    UsableElectromag(Super&&)      = delete;



    Super& super() { return *this; }
    Super const& super() const { return *this; }
    auto& operator*() { return super(); }
    auto& operator*() const { return super(); }


    UsableVecField<opts> E, B;
};




} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_ELECTROMAG_ELECTROMAG_FIXTURES_HPP */
