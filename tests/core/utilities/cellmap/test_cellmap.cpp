
#include <cstddef>
#include <iterator>
#include <string>
#include <vector>
#include <random>
#include <cmath>

#include "core/utilities/cellmap.hpp"
#include "core/utilities/bucketlist.hpp"
#include "core/utilities/box/box.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;

struct Obj
{
};

template<std::size_t dim>
struct Particle : BucketListItem
{
    std::array<int, dim> iCell;
    double delta;
    bool operator==(Particle const& other) const
    {
        bool equal = true;
        for (auto i = 0u; i < dim; ++i)
        {
            equal &= (iCell[i] == other.iCell[i]);
        }
        equal &= (std::abs(delta - other.delta) < 1e-12);
        return delta;
    }
};

TEST(CellMap, isDefaultConstructible)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    CellMap<dim, Obj, bucket_size> cm;
}




TEST(CellMap, sizeIsNbrOfElements)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 99u;
    CellMap<dim, Obj, bucket_size, int, Point<int, 2>> cm;
    Obj o1, o2;
    cm.addToCell(std::array<int, 2>{14, 27}, o1);
    EXPECT_EQ(cm.size({14, 27}), 1);
    cm.addToCell(Point{14, 26}, o2);
    EXPECT_EQ(cm.size({14, 26}), 1);
    EXPECT_EQ(cm.size({14, 25}), 0);
    EXPECT_EQ(cm.size({26, 14}), 0);
}


struct RemovableObj : BucketListItem
{
    RemovableObj(std::array<int, 2> cell)
        : iCell{cell}
    {
    }
    std::array<int, 2> iCell;
};

TEST(CellMap, itemCanBeRemoved)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 99u;
    CellMap<dim, RemovableObj, bucket_size, int, Point<int, 2>> cm;
    RemovableObj o1{{14, 27}}, o2{{14, 26}};
    cm.addToCell(std::array<int, 2>{14, 27}, o1);
    EXPECT_EQ(cm.size(), 1);
    cm.addToCell(Point{14, 26}, o2);
    EXPECT_EQ(cm.size(), 2);
    cm.erase(o2);
    EXPECT_EQ(1, cm.size());
}




TEST(CellMap, clearMakesSizeToZero)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 99u;
    CellMap<dim, Obj, bucket_size, int, Point<int, 2>> cm;
    Obj o1, o2;
    cm.addToCell(std::array<int, 2>{14, 27}, o1);
    cm.addToCell(Point{14, 26}, o2);
    cm.clear();
    EXPECT_EQ(0, cm.size());
}

template<std::size_t dim>
auto make_particles_in(PHARE::core::Box<int, dim> box, std::size_t nppc)
{
    std::vector<Particle<dim>> particles;
    particles.reserve(box.size() * nppc);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, 1.);
    for (auto const& cell : box)
    {
        for (auto ip = 0u; ip < nppc; ++ip)
        {
            Particle<dim> p;
            for (auto idim = 0u; idim < dim; ++idim)
            {
                p.iCell[idim] = cell[idim];
            }
            p.delta = dis(gen);
            particles.push_back(p);
        }
    }

    return particles;
}



TEST(CellMap, canStoreItemCollection)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    CellMap<dim, Particle<dim>, bucket_size> cm;
    Box<int, dim> patchbox{{10, 11}, {20, 21}};
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    EXPECT_EQ(patchbox.size() * nppc, particles.size());
    EXPECT_EQ(cm.size(), 0);
    EXPECT_EQ(cm.capacity(), 0);
    cm.add(particles);
    EXPECT_EQ(cm.size(), particles.size());
    EXPECT_EQ(cm.capacity(), cm.size());
}


TEST(CellMap, canStoreItemCollectionFromIterators)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    CellMap<dim, Particle<dim>, bucket_size> cm;
    Box<int, dim> patchbox{{10, 11}, {20, 21}};
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    EXPECT_EQ(patchbox.size() * nppc, particles.size());
    EXPECT_EQ(cm.size(), 0);
    EXPECT_EQ(cm.capacity(), 0);
    cm.add(std::begin(particles), std::end(particles));
    EXPECT_EQ(cm.size(), particles.size());
    EXPECT_EQ(cm.capacity(), cm.size());
}



TEST(CellMap, givesAccessToAllParticlesInACell)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, Particle<dim>, bucket_size, int, Point<int, 2>> cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    for (auto const& cell : patchbox)
    {
        auto& blist = cm.list_at(cell)->get();
        EXPECT_EQ(blist.size(), nppc);
        for (auto const& particle : blist)
        {
            EXPECT_EQ(Point{particle->iCell}, cell);
        }
    }
}



TEST(CellMap, emptyLeavesCapacityButZeroSize)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, Particle<dim>, bucket_size> cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    cm.empty();
    EXPECT_EQ(cm.size(), 0);
    EXPECT_TRUE(cm.is_empty());
    EXPECT_EQ(cm.capacity(), patchbox.size() * nppc);
}




TEST(CellMap, selectParticleInSubsetBox)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, Particle<dim>, bucket_size, int, Point<int, 2>> cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    Box<int, 2> selection_box{{12, 22}, {17, 23}};
    auto selected = cm.select(selection_box);
    for (auto const& particle : selected)
    {
        EXPECT_TRUE(isIn(Point{particle.iCell}, selection_box));
    }
    std::size_t isInCounter = 0;
    for (auto const& particle : particles)
    {
        if (isIn(Point{particle.iCell}, selection_box))
        {
            isInCounter++;
        }
    }
    EXPECT_EQ(isInCounter, selected.size());
}


TEST(CellMap, trimMemory)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, Particle<dim>, bucket_size, int, Point<int, dim>> cm;
    auto nppc      = 200u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    EXPECT_EQ(cm.size(), patchbox.size() * nppc);
    EXPECT_EQ(patchbox.size() * nppc, cm.capacity());

    cm.empty();
    EXPECT_EQ(cm.size(), 0);
    EXPECT_EQ(cm.capacity(), patchbox.size() * nppc);
    nppc      = 100;
    particles = make_particles_in(patchbox, nppc);
    cm.add(particles);
    EXPECT_EQ(cm.size(), patchbox.size() * nppc);
    EXPECT_EQ(patchbox.size() * nppc * 2, cm.capacity());
    EXPECT_EQ(0.5, cm.used_mem_ratio());
    cm.trim(1);
    EXPECT_EQ(patchbox.size() * nppc * 2, cm.capacity());
    cm.trim(0);
    EXPECT_EQ(patchbox.size() * nppc * 1, cm.capacity());
}


TEST(CellMap, emptyTrimsWithMaxEmpty1IfCapacityExceeds3TimesSize)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    using cellmap_t = CellMap<dim, Particle<dim>, bucket_size, int, Point<int, 2>>;
    cellmap_t cm;
    auto nppc      = 600u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);
    cm.empty();
    nppc      = 100u;
    particles = make_particles_in(patchbox, nppc);
    cm.add(particles);
    EXPECT_EQ(cm.capacity(), 600 * patchbox.size());
    cm.empty();
    EXPECT_EQ(cm.capacity(), particles.size() * 2);
}


TEST(CellMap, exportItems)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    using cellmap_t = CellMap<dim, Particle<dim>, bucket_size, int, Point<int, 2>>;
    cellmap_t cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);
    cm.empty();
    Box<int, 2> selectionBox{{14, 28}, {18, 37}};
    std::vector<Particle<dim>> selected;
    auto size = cm.size(selectionBox);
    selected.reserve(size);
    EXPECT_EQ(selected.size(), 0u);
    cm.export_to(selectionBox, selected);
    EXPECT_EQ(selected.size(), selected.capacity());
    for (auto const& p : selected)
    {
        EXPECT_TRUE(isIn(Point{p.iCell}, selectionBox));
    }
}

TEST(CellMap, threedim)
{
    auto constexpr dim         = 3u;
    auto constexpr bucket_size = 100;
    using cellmap_t            = CellMap<dim, Particle<dim>, bucket_size, int, Point<int, dim>>;
    Box<int, 3> patchbox{{10, 20, 30}, {25, 42, 54}};
    cellmap_t cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);
}


TEST(CellMap, trackParticle)
{
    auto constexpr dim         = 3u;
    auto constexpr bucket_size = 100;
    using cellmap_t            = CellMap<dim, Particle<dim>, bucket_size, int, Point<int, dim>>;
    Box<int, 3> patchbox{{10, 20, 30}, {25, 42, 54}};
    cellmap_t cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    EXPECT_EQ(cm.size(), particles.size());
    // pretends the particle change cell in x
    auto oldcell = particles[200].iCell;
    particles[200].iCell[0]++;

    cm.update(particles[200], oldcell);
    EXPECT_EQ(cm.size(), particles.size());

    for (auto const& cell : patchbox)
    {
        auto blist = cm.list_at(cell);
        if (blist)
        {
            for (auto const& part_ptr : blist->get())
            {
                EXPECT_TRUE(std::find(std::cbegin(particles), std::cend(particles), *part_ptr)
                            != std::cend(particles));
            }
        }
    }
}


template<std::size_t dim>
auto make_random_particles_in(Box<int, dim> const& box, std::size_t nppc)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::array<std::uniform_int_distribution<int>, dim> celldists;

    for (std::size_t idim = 0; idim < dim; ++idim)
        celldists[idim] = std::uniform_int_distribution<int>{box.lower[idim], box.upper[idim]};

    std::uniform_real_distribution<double> distdelta(0, 0.999999999999);

    std::vector<Particle<dim>> particles(nppc * box.size());
    for (std::size_t ip = 0; ip < particles.size(); ++ip)
    {
        for (std::size_t idim = 0; idim < dim; ++idim)
        {
            particles[ip].iCell[idim] = celldists[idim](gen);
        }
        particles[ip].delta = distdelta(gen);
    }
    return particles;
}


TEST(CellMap, sortArray)
{
    auto constexpr dim         = 3u;
    auto constexpr bucket_size = 100;
    using cellmap_t            = CellMap<dim, Particle<dim>, bucket_size, int, Point<int, dim>>;
    Box<int, 3> patchbox{{10, 20, 30}, {25, 42, 54}};
    cellmap_t cm;
    auto nppc      = 100u;
    auto particles = make_random_particles_in(patchbox, nppc);
    cm.add(particles);
    EXPECT_EQ(cm.size(), particles.size());

    std::size_t cell_jumps = 0;
    for (std::size_t ipart = 1; ipart < particles.size(); ++ipart)
    {
        if (particles[ipart].iCell > particles[ipart - 1].iCell)
        {
            cell_jumps++;
        }
    }
    EXPECT_NE(patchbox.size() - 1, cell_jumps);
    std::vector<Particle<3>> sorted(particles.size());
    get_sorted(cm, patchbox, sorted);

    cell_jumps = 0;
    for (std::size_t ipart = 1; ipart < sorted.size(); ++ipart)
    {
        if (sorted[ipart].iCell > sorted[ipart - 1].iCell)
        {
            cell_jumps++;
        }
    }
    EXPECT_EQ(patchbox.size() - 1, cell_jumps);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

