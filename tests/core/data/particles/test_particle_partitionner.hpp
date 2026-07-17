// tests/core/data/particles/test_particle_partitionner.cpp

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_utilities.hpp"
#include "core/data/particles/particle_array_partitionner.hpp"
#include "core/utilities/point/point.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "tests/core/data/particles/test_particles.hpp"


using namespace PHARE::core;

class ParticlePartitioner : public ::testing::Test
{
protected:
    Particle<3> part;

public:
    ParticlePartitioner()
        : part{0.01, 1., {43, 75, 92}, {0.002, 0.2, 0.8}, {1.8, 1.83, 2.28}}
    {
    }
};

TEST_F(ParticlePartitioner, box_remove_1)
{
    using box_t = Box<int, 1>;
    box_t a{{0}, {2}};
    box_t b{{1}, {1}};
    auto const remaining = a.remove(b);
    assert(not any_overlaps(remaining));
    EXPECT_EQ(2, sum_from(remaining, [](auto& r) { return r.size(); }));
}

TEST_F(ParticlePartitioner, box_remove_2)
{
    using box_t = Box<int, 2>;
    box_t a{{0, 0}, {2, 2}};
    box_t b{{1, 1}, {1, 1}};
    auto const remaining = a.remove(b);
    assert(not any_overlaps(remaining));
    EXPECT_EQ(std::pow(3, 2) - 1, sum_from(remaining, [](auto& r) { return r.size(); }));
}


TEST_F(ParticlePartitioner, box_remove_3)
{
    auto constexpr static dim = 3;
    using box_t               = Box<int, dim>;
    using point_t             = Point<int, dim>;

    {
        box_t a{{0, 0, 0}, {4, 4, 4}};
        box_t b{{0, 1, 0}, {4, 4, 4}};
        box_t r{{0, 0, 0}, {4, 0, 4}};
        auto const remaining = a.remove(b);
        assert(not any_overlaps(remaining));
        assert(remaining[0] == r);
    }

    auto expect_cells = [](auto const& boxes, auto const& skip_point) {
        for (int x = 0; x < 3; ++x)
            for (int y = 0; y < 3; ++y)
                for (int z = 0; z < 3; ++z)
                    if (point_t p{x, y, z}; p != skip_point)
                        assert(isIn(p, boxes));
    };

    {
        point_t p{0, 0, 0};
        box_t a{{0, 0, 0}, {2, 2, 2}};
        auto const remaining = a.remove(box_t{p, p});
        assert(not any_overlaps(remaining));
        EXPECT_EQ(std::pow(3, 3) - 1, sum_from(remaining, [](auto& r) { return r.size(); }));
        expect_cells(remaining, p);
    }

    {
        point_t p{0, 0, 2};
        box_t a{{0, 0, 0}, {2, 2, 2}};
        auto const remaining = a.remove(box_t{p, p});
        assert(not any_overlaps(remaining));
        EXPECT_EQ(std::pow(3, 3) - 1, sum_from(remaining, [](auto& r) { return r.size(); }));
        expect_cells(remaining, p);
    }

    {
        point_t p{2, 2, 2};
        box_t a{{0, 0, 0}, {2, 2, 2}};
        auto const remaining = a.remove(box_t{p, p});
        assert(not any_overlaps(remaining));
        EXPECT_EQ(std::pow(3, 3) - 1, sum_from(remaining, [](auto& r) { return r.size(); }));
        expect_cells(remaining, p);
    }

    {
        point_t p{1, 1, 1};
        box_t a{{0, 0, 0}, {2, 2, 2}};
        auto const remaining = a.remove(box_t{p, p});
        assert(not any_overlaps(remaining));
        EXPECT_EQ(std::pow(3, 3) - 1, sum_from(remaining, [](auto& r) { return r.size(); }));
        expect_cells(remaining, p);
    }
}



TEST_F(ParticlePartitioner, simple_unique_overlaps_9)
{
    auto constexpr static dim               = 3;
    auto constexpr static extra_ghost_cells = 1;
    using box_t                             = Box<int, dim>;

    box_t middle_box{{3, 3, 3}, {5, 5, 5}};

    std::vector<box_t> neighbor_boxes{
        // z == 0
        {{3, 3, 0}, {5, 5, 2}},
        // // z = 3
        {{3, 6, 3}, {5, 8, 5}},
    };
    assert(not any_overlaps(neighbor_boxes));
    assert(not any_overlaps(neighbor_boxes, middle_box));

    auto neighbor_ghost_boxes
        = generate_from([](auto box) { return box.grow(extra_ghost_cells); }, neighbor_boxes);
    assert(all_overlaps(neighbor_ghost_boxes, middle_box));

    auto overlaps = distinct_overlaps(neighbor_ghost_boxes, middle_box);
    assert(not any_overlaps(overlaps));

    EXPECT_EQ(15, sum_from(overlaps, [](auto& r) { return r.size(); }));
}

TEST_F(ParticlePartitioner, simple_unique_overlaps)
{
    auto constexpr static dim               = 3;
    auto constexpr static extra_ghost_cells = 1;
    using box_t                             = Box<int, dim>;

    box_t middle_box{{3, 3, 3}, {5, 5, 5}};

    std::vector<box_t> neighbor_boxes{
        // z == 0
        {{3, 3, 0}, {5, 5, 2}},
        // // z = 3
        {{3, 0, 3}, {5, 2, 5}},
        {{0, 3, 3}, {2, 5, 5}},
        {{6, 3, 3}, {8, 5, 5}},
        {{3, 6, 3}, {5, 8, 5}},
        // // z = 6
        {{3, 3, 6}, {5, 5, 8}},
    };
    assert(not any_overlaps(neighbor_boxes));
    assert(not any_overlaps(neighbor_boxes, middle_box));

    auto neighbor_ghost_boxes
        = generate_from([](auto box) { return box.grow(extra_ghost_cells); }, neighbor_boxes);
    assert(all_overlaps(neighbor_ghost_boxes, middle_box));

    auto overlaps = distinct_overlaps(neighbor_ghost_boxes, middle_box);
    assert(not any_overlaps(overlaps));

    EXPECT_EQ(std::pow(3, 3) - 1, sum_from(overlaps, [](auto& r) { return r.size(); }));
}


TEST_F(ParticlePartitioner, full_overlaps)
{
    auto constexpr static dim               = 3;
    auto constexpr static extra_ghost_cells = 1;
    using box_t                             = Box<int, dim>;

    box_t middle_box{{3, 3, 3}, {5, 5, 5}};

    std::vector<box_t> neighbor_boxes{
        // z == 0
        {{0, 0, 0}, {2, 2, 2}},
        {{3, 0, 0}, {5, 2, 2}},
        {{6, 0, 0}, {8, 2, 2}},

        {{0, 3, 0}, {2, 5, 2}},
        {{3, 3, 0}, {5, 5, 2}},
        {{6, 3, 0}, {8, 5, 2}},

        {{0, 6, 0}, {2, 8, 2}},
        {{3, 6, 0}, {5, 8, 2}},
        {{6, 6, 0}, {8, 8, 2}},

        // // z = 3
        {{0, 0, 3}, {2, 2, 5}},
        {{3, 0, 3}, {5, 2, 5}},
        {{6, 0, 3}, {8, 2, 5}},

        {{0, 3, 3}, {2, 5, 5}},
        // {{3, 3, 3}, {5, 5, 5}}, // skip
        {{6, 3, 3}, {8, 5, 5}},

        {{0, 6, 3}, {2, 8, 5}},
        {{3, 6, 3}, {5, 8, 5}},
        {{6, 6, 3}, {8, 8, 5}},

        // // z = 6
        {{0, 0, 6}, {2, 2, 8}},
        {{3, 0, 6}, {5, 2, 8}},
        {{6, 0, 6}, {8, 2, 8}},

        {{0, 3, 6}, {2, 5, 8}},
        {{3, 3, 6}, {5, 5, 8}},
        {{6, 3, 6}, {8, 5, 8}},

        {{0, 6, 6}, {2, 8, 8}},
        {{3, 6, 6}, {5, 8, 8}},
        {{6, 6, 6}, {8, 8, 8}},
    };

    assert(not any_overlaps(neighbor_boxes));
    assert(not any_overlaps(neighbor_boxes, middle_box));

    auto neighbor_ghost_boxes
        = generate_from([](auto box) { return box.grow(extra_ghost_cells); }, neighbor_boxes);

    assert(all_overlaps(neighbor_ghost_boxes, middle_box));

    auto overlaps = distinct_overlaps(neighbor_ghost_boxes, middle_box);
    assert(not any_overlaps(overlaps));

    EXPECT_EQ(std::pow(3, 3) - 1, sum_from(overlaps, [](auto& r) { return r.size(); }));
}


TEST_F(ParticlePartitioner, full_overlaps_5_x_5)
{
    auto constexpr static dim               = 3;
    auto constexpr static extra_ghost_cells = 2;
    using box_t                             = Box<int, dim>;
    using point_t                           = Point<int, dim>;

    box_t middle_box{{10, 10, 10}, {14, 14, 14}};

    std::vector<box_t> neighbor_boxes;
    for (int x = 1; x < 4; ++x)
    {
        auto x0 = x * 5;

        for (int y = 1; y < 4; ++y)
        {
            auto y0 = y * 5;

            for (int z = 1; z < 4; ++z)
            {
                auto z0 = z * 5;

                point_t p{x0, y0, z0};
                if (p == middle_box.lower)
                    continue;

                neighbor_boxes.push_back(box_t{{x0, y0, z0}, {x0 + 4, y0 + 4, z0 + 4}});
            }
        }
    }

    assert(not any_overlaps(neighbor_boxes));
    assert(not any_overlaps(neighbor_boxes, middle_box));

    auto neighbor_ghost_boxes
        = generate_from([](auto box) { return box.grow(extra_ghost_cells); }, neighbor_boxes);

    assert(all_overlaps(neighbor_ghost_boxes, middle_box));

    auto overlaps = distinct_overlaps(neighbor_ghost_boxes, middle_box);
    assert(not any_overlaps(overlaps));

    EXPECT_EQ(std::pow(5, 3) - 1, sum_from(overlaps, [](auto& r) { return r.size(); }));
}


TEST_F(ParticlePartitioner, partition_ghosts)
{
    // here we are only separating particles into three
    // domain
    // ghost layer
    // outside ghost layer

    auto constexpr static dim               = 3;
    auto constexpr static extra_ghost_cells = 2;
    using box_t                             = Box<int, dim>;

    auto middle_box       = box_t{{10, 10, 10}, {14, 14, 14}};
    auto super_ghost_box  = grow(middle_box, extra_ghost_cells + 1);
    auto middle_ghost_box = grow(middle_box, extra_ghost_cells);
    auto ghost_boxes      = middle_ghost_box.remove(middle_box);

    auto L = [&](auto&& particles) {
        using ParticleArray_t = std::decay_t<decltype(particles)>;
        using Partitioner = ParticleArrayPartitioner<ParticleArray_t::alloc_mode, ParticleArray_t>;

        std::size_t ppc = 10;
        add_particles_in(particles, super_ghost_box, ppc);
        assert(particles.size() == std::pow(11, 3) * ppc);

        auto iterators
            = Partitioner{particles}(std::array{middle_box, grow(middle_box, extra_ghost_cells)});

        assert(iterators.size() == 2);
        assert(particles.begin() == iterators[0].begin());
        assert(std::distance(particles.begin(), iterators[0].end()) == std::pow(5, 3) * ppc);
        assert(iterators[0].size() == std::pow(5, 3) * ppc);

        for (auto p : iterators[0])
        {
            assert((not isIn(p.iCell(), ghost_boxes)));
        }

        assert(iterators[1].size() == (std::pow(9, 3) - std::pow(5, 3)) * ppc);
        for (auto it = iterators[1].begin(); it != iterators[1].end(); ++it)
        {
            assert(isIn((*it).iCell(), ghost_boxes));
        }

        auto outside_ghost_box = std::distance(iterators.back().end(), particles.end());
        assert(outside_ghost_box == (std::pow(11, 3) - std::pow(9, 3)) * ppc);
        for (auto it = iterators[1].end(); it != particles.end(); ++it)
        {
            assert(not isIn((*it).iCell(), ghost_boxes) and isIn((*it).iCell(), super_ghost_box));
        }
    };

    L(AoSParticleArray<dim>{});
    PHARE_WITH_THRUST(L(SoAParticleArray<dim>{}));
}



// TEST_F(ParticlePartitioner, partition_overlaps)
// {
//     auto constexpr static dim               = 3;
//     auto constexpr static extra_ghost_cells = 2;
//     using box_t                             = Box<int, dim>;
//     using point_t                           = Point<int, dim>;

//     box_t middle_box{{10, 10, 10}, {14, 14, 14}};

//     std::vector<box_t> neighbor_boxes;
//     for (int x = 1; x < 4; ++x)
//     {
//         auto x0 = x * 5;

//         for (int y = 1; y < 4; ++y)
//         {
//             auto y0 = y * 5;

//             for (int z = 1; z < 4; ++z)
//             {
//                 auto z0 = z * 5;

//                 point_t p{x0, y0, z0};
//                 if (p == middle_box.lower)
//                     continue;

//                 neighbor_boxes.push_back(box_t{{x0, y0, z0}, {x0 + 4, y0 + 4, z0 + 4}});
//             }
//         }
//     }

//     assert(not any_overlaps(neighbor_boxes));
//     assert(not any_overlaps(neighbor_boxes, middle_box));
//     auto neighbor_ghost_boxes
//         = generate_from([](auto box) { return box.grow(extra_ghost_cells); }, neighbor_boxes);
//     assert(all_overlaps(neighbor_ghost_boxes, middle_box));

//     auto overlaps = distinct_overlaps(neighbor_ghost_boxes, middle_box);
//     assert(not any_overlaps(overlaps));
//     EXPECT_EQ(std::pow(5, 3) - 1, sum_from(overlaps, [](auto& r) { return r.size(); }));

//     auto middle_ghost_box = grow(middle_box, extra_ghost_cells);
//     auto ghost_boxes      = middle_ghost_box.remove(middle_box);

//     auto L = [&](auto&& particles) {
//         std::size_t ppc = 10;
//         add_particles_in(particles, middle_ghost_box, ppc);
//         assert(particles.size() == std::pow(9, 3) * ppc);

//         auto all_boxes = neighbor_boxes;
//         neighbor_boxes.emplace_back();
//         // auto iterators = partition<extra_ghost_cells>(particles, middle_box, neighbor_boxes);
//         auto iterators = ParticleArrayPartitioner<ParticleArray_t>{particles}(
//             std::array{middle_box, grow(middle_box, extra_ghost_cells)});

//         assert(iterators.size() > 1);

//         assert(particles.begin() == iterators[0].begin());
//         assert(std::distance(particles.begin(), iterators[0].end()) == ppc);

//         for (auto it = particles.begin(); it != iterators[0].end(); ++it)
//             assert((not isIn((*it).iCell(), ghost_boxes)) and (not isIn((*it).iCell(),
//             overlaps)));


//         for (std::size_t i = 1; i < iterators.size(); ++i)
//             for (auto it = iterators[i].begin(); it != iterators[i].end(); ++it)
//                 assert(isIn((*it).iCell(), overlaps));

//         for (auto it = iterators.back().end(); it != particles.end(); ++it)
//             assert(isIn((*it).iCell(), ghost_boxes));

//         auto n_ghost_particles = std::distance(iterators.back().end(), particles.end());
//         assert(n_ghost_particles == (std::pow(9, 3) - std::pow(5, 3)) * ppc);
//     };

//     L(AoSParticleArray<dim>{});
//     // L(SoAParticleArray<dim>{});
// }




// int main(int argc, char** argv)
//{
//     ::testing::InitGoogleTest(&argc, argv);
//
//     return RUN_ALL_TESTS();
// }
