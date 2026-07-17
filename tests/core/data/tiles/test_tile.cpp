
#include "test_tile.hpp"

#include <unordered_set>


TEST(TraverseTiles, visitsEachTileExactlyOnce)
{
    Box<int, 3> const box{{0, 0, 0}, {23, 23, 23}};
    TileSet<TileMock<3>> ts{box};
    TileSet<TileMock<3>>::build_links(ts);

    std::unordered_set<TileMock<3>*> visited;
    traverse_tiles(ts, box, [&](auto& tile) {
        EXPECT_EQ(visited.count(&tile), 0u);
        visited.insert(&tile);
    });
    EXPECT_EQ(visited.size(), ts.size());
}


TEST(TraverseTiles, fuzzVariousShapesAndSubBoxes)
{
    // box upper per dim: 12 → uniform [6,6], 16 → [6,6,4], 20 → [6,6,4,4]
    for (int ux : {12, 16, 20})
        for (int uy : {12, 16, 20})
            for (int uz : {12, 16, 20})
            {
                Box<int, 3> const box{{0, 0, 0}, {ux - 1, uy - 1, uz - 1}};
                TileSet<TileMock<3>> ts{box};
                TileSet<TileMock<3>>::build_links(ts);

                for (int sx = 2; sx <= 10; ++sx)
                    for (int sy = 2; sy <= 10; ++sy)
                        for (int sz = 2; sz <= 10; ++sz)
                        {
                            Box<int, 3> const sub_box{{0, 0, 0}, {sx, sy, sz}};
                            SCOPED_TRACE("box=" + std::to_string(ux) + "x" + std::to_string(uy)
                                         + "x" + std::to_string(uz)
                                         + " sub=" + std::to_string(sx + 1) + "x"
                                         + std::to_string(sy + 1) + "x" + std::to_string(sz + 1));

                            std::unordered_set<TileMock<3>*> visited;
                            traverse_tiles(ts, sub_box, [&](auto& tile) {
                                EXPECT_EQ(visited.count(&tile), 0u);
                                visited.insert(&tile);
                            });

                            for (auto& tile : ts)
                            {
                                bool const overlaps = (sub_box * tile).has_value();
                                EXPECT_EQ(visited.count(&tile), overlaps ? 1u : 0u);
                            }

                            EXPECT_GT(visited.size(), 0);
                        }
            }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
