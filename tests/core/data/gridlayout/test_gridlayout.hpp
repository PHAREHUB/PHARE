#ifndef TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/utilities/types.hpp"


template<typename GridLayout>
class TestGridLayout : public GridLayout
{ // to expose a default constructor
public:
    auto static constexpr dim = GridLayout::dimension;

    TestGridLayout() = default;

    TestGridLayout(std::uint32_t const cells)
        : GridLayout{PHARE::core::ConstArray<double, dim>(1.0 / cells),
                     PHARE::core::ConstArray<std::uint32_t, dim>(cells),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}}
    {
    }


    TestGridLayout(std::array<double, dim> const dl, auto&&... args)
        : GridLayout{dl, args...}
    {
    }

    auto static make(std::uint32_t cells) { return TestGridLayout{cells}; }
    auto static make(PHARE::core::Box<int, dim> const box)
    {
        using namespace PHARE::core;
        auto const shape = for_N_make_array<dim>(
            [&](auto i) -> std::uint32_t { return box.upper[i] - box.lower[i] + 1; });
        auto const dl = for_N_make_array<dim>([&](auto i) -> double { return 1.0 / shape[i]; });
        auto const origin
            = for_N_make_array<dim>([&](auto i) -> double { return box.lower[i] * dl[i]; });

        return TestGridLayout{dl, shape, origin, box};
    }
};



#endif /*TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP*/
