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

    TestGridLayout(std::uint32_t cells)
        : GridLayout{PHARE::core::ConstArray<double, dim>(1.0 / cells),
                     PHARE::core::ConstArray<std::uint32_t, dim>(cells),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}}
    {
    }

    auto static make(std::uint32_t cells) { return TestGridLayout{cells}; }
};

#endif /*TESTS_CORE_DATA_GRIDLAYOUT_TEST_GRIDLAYOUT_HPP*/
