#ifndef TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP

#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include <cstdint>
#include <type_traits>


template<typename GridLayout>
class TestGridLayout : public GridLayout
{ // to expose a default constructor
public:
    auto static constexpr dim = GridLayout::dimension;

    TestGridLayout() = default;

    TestGridLayout(std::uint32_t const& cells)
        : GridLayout{PHARE::core::ConstArray<double, dim>(1.0 / cells),
                     PHARE::core::ConstArray<std::uint32_t, dim>(cells),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}}
    {
    }


    TestGridLayout(PHARE::core::Box<int, dim> const& amrbox)
        : GridLayout{PHARE::core::ConstArray<double, dim>(.1),
                     amrbox.shape().template toArray<std::uint32_t>(),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)},
                     amrbox}
    {
    }

    template<typename... Args>
    TestGridLayout(Args&&... args)
        requires(std::is_constructible_v<GridLayout, Args...>)
        : GridLayout{args...}
    {
    }

    auto static make(auto const& cells) { return TestGridLayout{cells}; }
};

#endif /*TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP*/
