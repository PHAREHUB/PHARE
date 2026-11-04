// Sorted mode of ParticleArray: bag sort (sortBags), row layout (layoutRows), compaction and the
// box queries / cell counter on sorted arrays, checked against brute force.

#include "core/data/particles/particle_array.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <cstddef>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>


namespace PHARE::core
{
template<std::size_t dim>
bool inBox(std::array<int, dim> const& cell, Box<int, dim> const& box)
{
    for (std::size_t d = 0; d < dim; ++d)
        if (cell[d] < box.lower[d] or cell[d] > box.upper[d])
            return false;
    return true;
}


/** random "pushes" on a bag layout: 5 % of the live particles move to a neighbour cell (some
 * leave the owned box), particles are appended (sometimes all in one cell: row overflow and
 * cascades), then sortBags. After each sort: rows hold their key only, the live particles are
 * the expected ones (ids in the weight), the box queries match brute force; every 37 steps a
 * compaction must leave a packed CSR of the same particles.
 */
template<std::size_t dim>
void randomBagSorts(unsigned const seed, int const steps)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> uni(0, 1 << 30);

    Box<int, dim> box;
    for (std::size_t d = 0; d < dim; ++d)
    {
        box.lower[d] = -2;
        box.upper[d] = 5 + 2 * static_cast<int>(d);
    }
    Box<int, dim> owned = box; // a ring of one non owned cell
    for (std::size_t d = 0; d < dim; ++d)
    {
        owned.lower[d] += 1;
        owned.upper[d] -= 1;
    }

    ParticleArray<dim> array(box);
    int const K = array.nbrCells();
    auto keyOf  = [&](Particle<dim> const& p) {
        return inBox<dim>(p.iCell, owned) ? array.key(p.iCell) : K;
    };

    double nextId = 1;
    auto make     = [&]() {
        Particle<dim> p;
        p.weight = nextId++;
        for (std::size_t d = 0; d < dim; ++d)
            p.iCell[d] = owned.lower[d] + uni(gen) % (owned.upper[d] - owned.lower[d] + 1);
        return p;
    };

    std::multiset<double> expected;
    for (int i = 0; i < 20 * K; ++i)
    {
        auto const p = make();
        array.push_back(p);
        expected.insert(p.weight);
    }

    typename ParticleArray<dim>::SortScratch scratch;
    std::vector<std::pair<std::size_t, int>> crossers;

    Box<int, dim> query;
    for (std::size_t d = 0; d < dim; ++d)
    {
        query.lower[d] = 0;
        query.upper[d] = 2 + static_cast<int>(d);
    }

    auto verify = [&](int const step) {
        SCOPED_TRACE("dim " + std::to_string(dim) + " seed " + std::to_string(seed) + " step "
                     + std::to_string(step));
        ASSERT_TRUE(array.checkBags(keyOf));
        ASSERT_TRUE(array.is_consistent());

        std::multiset<double> live;
        std::size_t inQuery = 0;
        for (int k = 0; k < K; ++k)
            for (int i = array.csr()[k]; i < array.rowEnds()[k]; ++i)
            {
                live.insert(array[i].weight);
                inQuery += inBox<dim>(array[i].iCell, query);
            }
        EXPECT_EQ(live, expected);
        EXPECT_EQ(array.liveSize(), expected.size());

        std::size_t runs = 0;
        array.forEachLiveRun([&](std::size_t first, std::size_t last) { runs += last - first; });
        EXPECT_EQ(runs, expected.size());

        ParticleArray<dim> exported(box);
        array.export_particles(query, exported);
        auto const counter   = array.cellCounter(query);
        std::size_t perCells = 0;
        for_each_point(query, [&](auto const& cell) { perCells += counter(cell.toArray()); });
        EXPECT_EQ(array.nbr_particles_in(query), inQuery);
        EXPECT_EQ(exported.size(), inQuery);
        EXPECT_EQ(perCells, inQuery);
    };

    array.sortBags(keyOf, 0, crossers, scratch);
    verify(-1);

    for (int step = 0; step < steps; ++step)
    {
        if (!array.bagged()) // after a compaction: lay out again
            array.sortBags(keyOf, 0, crossers, scratch);

        crossers.clear();
        std::size_t const oldPrefix = array.csrPrefix();
        for (int k = 0; k < K; ++k)
            for (int i = array.csr()[k]; i < array.rowEnds()[k]; ++i)
                if (uni(gen) % 100 < 5)
                {
                    auto& p        = array[i];
                    auto const d   = static_cast<std::size_t>(uni(gen)) % dim;
                    int const sign = uni(gen) % 2 ? 1 : -1;
                    p.iCell[d] += sign;
                    if (!inBox<dim>(p.iCell, box))
                    {
                        p.iCell[d] -= sign;
                        continue;
                    }
                    crossers.emplace_back(i, k);
                    if (keyOf(p) == K)
                        expected.erase(expected.find(p.weight));
                }

        int const nAppend = uni(gen) % 60;
        for (int i = 0; i < nAppend; ++i)
        {
            auto p = make();
            if (step % 7 == 0)
                p.iCell = owned.lower.toArray();
            array.push_back(p);
            expected.insert(p.weight);
        }

        array.sortBags(keyOf, oldPrefix, crossers, scratch);
        for (auto const& p : scratch.leavers)
            EXPECT_EQ(keyOf(p), K);
        verify(step);

        if (step % 37 == 36)
        {
            array.compact();
            EXPECT_FALSE(array.bagged());
            EXPECT_TRUE(array.checkCSR(keyOf, array.size()));
            std::multiset<double> live;
            for (auto const& p : array)
                live.insert(p.weight);
            EXPECT_EQ(live, expected);
        }
    }
}

TEST(ParticleArraySorted, randomBagSorts1D)
{
    for (unsigned seed = 0; seed < 4; ++seed)
        randomBagSorts<1>(seed, 300);
}
TEST(ParticleArraySorted, randomBagSorts2D)
{
    for (unsigned seed = 0; seed < 4; ++seed)
        randomBagSorts<2>(seed, 300);
}
TEST(ParticleArraySorted, randomBagSorts3D)
{
    for (unsigned seed = 0; seed < 4; ++seed)
        randomBagSorts<3>(seed, 150);
}



/** CellCounter == nbr_particles_in(cell) on an unsorted (cell map) array, on a row layout whose
 * unsorted appends include cells outside the array box, and with more appends after the layout;
 * the query box sticks out of the array box.
 */
template<std::size_t dim>
void cellCounterMatchesCellQueries(unsigned const seed)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> coord(-2, 9);

    Box<int, dim> box;
    for (std::size_t d = 0; d < dim; ++d)
    {
        box.lower[d] = 0;
        box.upper[d] = 7 + static_cast<int>(d);
    }
    ParticleArray<dim> array(box);
    auto make = [&](bool const inside) {
        Particle<dim> p;
        for (std::size_t d = 0; d < dim; ++d)
        {
            int const c = coord(gen);
            p.iCell[d]  = inside ? std::clamp(c, 0, 7) : c;
        }
        return p;
    };
    for (int i = 0; i < 3000; ++i)
        array.push_back(make(true));

    auto check = [&](std::string const& what) {
        SCOPED_TRACE("dim " + std::to_string(dim) + " seed " + std::to_string(seed) + " " + what);
        Box<int, dim> query = box;
        Box<int, dim> sweep = box;
        for (std::size_t d = 0; d < dim; ++d)
            query.lower[d] -= 1;
        if (array.csrValid()) // the cell map only knows the cells of the box
            for (std::size_t d = 0; d < dim; ++d)
            {
                sweep.lower[d] -= 2;
                sweep.upper[d] += 2;
            }
        auto const counter = array.cellCounter(query);
        for_each_point(sweep, [&](auto const& cell) {
            EXPECT_EQ(counter(cell.toArray()), array.nbr_particles_in(cell.toArray()));
        });
    };

    check("cell map");

    typename ParticleArray<dim>::SortScratch scratch;
    int const K = array.nbrCells();
    array.layoutRows(
        [&](auto const& p) { return inBox<dim>(p.iCell, box) ? array.key(p.iCell) : K; }, scratch);
    ASSERT_TRUE(array.csrValid());
    check("rows");

    for (int i = 0; i < 500; ++i)
        array.push_back(make(false));
    check("rows + appends");
}

TEST(ParticleArraySorted, cellCounterMatchesCellQueries)
{
    for (unsigned seed = 0; seed < 3; ++seed)
    {
        cellCounterMatchesCellQueries<1>(seed);
        cellCounterMatchesCellQueries<2>(seed);
        cellCounterMatchesCellQueries<3>(seed);
    }
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
