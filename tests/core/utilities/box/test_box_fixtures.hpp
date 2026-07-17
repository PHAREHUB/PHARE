#ifndef PHARE_TEST_CORE_UTILITIES_BOX_TEST_BOX_FIXTURES_HPP
#define PHARE_TEST_CORE_UTILITIES_BOX_TEST_BOX_FIXTURES_HPP

#include "phare_core.hpp"
#include "core/utilities/box/box.hpp"

#include <array>
#include <tuple>
#include <vector>

namespace PHARE::core
{



// get middle box and 8 adjacent neighbour boxes
template<std::size_t extra_ghost_cells_ = 2>
auto get_some_boxes(std::size_t const size = 5)
{
    abort_if(size < 3); // no

    using T                                 = int;
    std::size_t constexpr static dim        = 3; // this is only 3d
    auto constexpr static extra_ghost_cells = extra_ghost_cells_;
    using box_t                             = Box<T, dim>;
    using point_t                           = Point<T, dim>;

    T sn1 = size - 1;
    T l   = size;
    T u   = size + sn1;
    box_t middle_box{{l, l, l}, {u, u, u}};

    std::vector<box_t> neighbor_boxes;
    for (T x = 0; x < 3; ++x)
    {
        T x0 = x * size;
        for (T y = 0; y < 3; ++y)
        {
            T y0 = y * size;
            for (T z = 0; z < 3; ++z)
            {
                T z0 = z * size;
                if (point_t p{x0, y0, z0}; p == middle_box.lower)
                    continue;
                neighbor_boxes.push_back(box_t{{x0, y0, z0}, {x0 + sn1, y0 + sn1, z0 + sn1}});
            }
        }
    }

    assert(neighbor_boxes.size() == 26);

    return std::make_tuple(middle_box, neighbor_boxes);
}

} // namespace PHARE::core

#endif /*PHARE_TEST_CORE_UTILITIES_BOX_TEST_BOX_FIXTURES_HPP*/