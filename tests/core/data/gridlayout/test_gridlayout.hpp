
#ifndef TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"


template<typename GridLayout>
class TestGridLayout : public GridLayout // to expose a default constructor
{
    auto constexpr static dimension    = GridLayout::dimension;
    double constexpr static level_0_dx = .1;

public:
    using Super = GridLayout;

    auto static constexpr dim = GridLayout::dimension;

    TestGridLayout() = default;

    TestGridLayout(double const dl, std::uint32_t const cells)
        : GridLayout{PHARE::core::ConstArray<double, dim>(dl),
                     PHARE::core::ConstArray<std::uint32_t, dim>(cells),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}}
    {
    }

    TestGridLayout(std::uint32_t const cells)
        : TestGridLayout{.1, cells}
    {
    }


    TestGridLayout(std::array<double, dim> const dl, auto&&... args)
        : GridLayout{dl, args...}
    {
    }

    TestGridLayout(PHARE::core::Box<int, dim> const& amrbox, double const dl = .1)
        : GridLayout{PHARE::core::ConstArray<double, dim>(dl),
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

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }
};



#endif /*TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP*/