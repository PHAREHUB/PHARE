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
    using Super               = GridLayout;
    auto static constexpr dim = GridLayout::dimension;

    TestGridLayout() = default;

    TestGridLayout(double dl, std::uint32_t cells)
        : GridLayout{PHARE::core::ConstArray<double, dim>(dl),
                     PHARE::core::ConstArray<std::uint32_t, dim>(cells),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}}
    {
    }

    TestGridLayout(std::uint32_t const cells)
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

    auto static make(std::uint32_t cells) { return TestGridLayout{cells}; }
    auto static make(double dl, std::uint32_t cells) { return TestGridLayout{dl, cells}; }
    auto static make(PHARE::core::Box<int, dim> const& amrbox) { return TestGridLayout{amrbox}; }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }
};

#endif /*TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP*/
