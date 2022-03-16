
#include <cstddef>
#include <iterator>
#include <string>
#include <vector>
#include <random>
#include <cmath>

#include "core/utilities/cellmap.hpp"
#include "core/utilities/bucketlist.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;


struct Obj
{
};

template<std::size_t dim>
struct Particle
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
    CellMap<dim, bucket_size> cm;
}




TEST(CellMap, sizeIsNbrOfElements)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 99u;
    CellMap<dim, bucket_size, int, Point<int, 2>> cm;
    cm.addToCell(std::array<int, 2>{14, 27}, 0);
    EXPECT_EQ(cm.size({14, 27}), 1);
    cm.addToCell(Point{14, 26}, 1);
    EXPECT_EQ(cm.size({14, 26}), 1);
    EXPECT_EQ(cm.size({14, 25}), 0);
    EXPECT_EQ(cm.size({26, 14}), 0);
}

TEST(CellMap, canSortCells)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 99u;
    CellMap<dim, bucket_size, int, Point<int, 2>> cm;

    // add indexes not in order
    // in two different cells
    cm.addToCell(std::array<int, 2>{14, 27}, 1);
    cm.addToCell(Point{14, 27}, 0);
    cm.addToCell(Point{14, 27}, 3);
    cm.addToCell(Point{14, 27}, 2);

    cm.addToCell(std::array<int, 2>{12, 22}, 3);
    cm.addToCell(Point{12, 22}, 2);
    cm.addToCell(Point{12, 22}, 0);
    cm.addToCell(Point{12, 22}, 1);

    cm.sort();

    for (auto const& [cell, indexes] : cm)
    {
        int i = -1;
        for (auto index : indexes)
        {
            if (i == -1)
            {
                i = index;
            }
            else
            {
                EXPECT_GT(index, i);
                i = index;
            }
        }
    }
}


TEST(CellMap, itemCanBeRemoved)
{
    std::array<Particle<2>, 2> particles;
    particles[0].iCell[0]      = 14;
    particles[0].iCell[1]      = 27;
    particles[1].iCell[0]      = 14;
    particles[1].iCell[1]      = 26;
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 99u;
    CellMap<dim, bucket_size, int, Point<int, 2>> cm;
    cm.addToCell(std::array<int, 2>{14, 27}, 0);
    EXPECT_EQ(cm.size(), 1);
    cm.addToCell(Point{14, 26}, 1);
    EXPECT_EQ(cm.size(), 2);
    cm.erase(particles, 1);
    EXPECT_EQ(1, cm.size());
}




TEST(CellMap, clearMakesSizeToZero)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 99u;
    CellMap<dim, bucket_size, int, Point<int, 2>> cm;
    cm.addToCell(std::array<int, 2>{14, 27}, 0);
    cm.addToCell(Point{14, 26}, 1);
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
    CellMap<dim, bucket_size> cm;
    Box<int, dim> patchbox{{10, 11}, {20, 21}};
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    EXPECT_EQ(patchbox.size() * nppc, particles.size());
    EXPECT_EQ(cm.size(), 0);
    EXPECT_EQ(cm.capacity(), 0);
    cm.add(particles);
    EXPECT_EQ(cm.size(), particles.size());
    // EXPECT_EQ(cm.capacity(), cm.size()); // not true with default vector
}


TEST(CellMap, canStoreItemCollectionFromIterators)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    CellMap<dim, bucket_size> cm;
    Box<int, dim> patchbox{{10, 11}, {20, 21}};
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    EXPECT_EQ(patchbox.size() * nppc, particles.size());
    EXPECT_EQ(cm.size(), 0);
    EXPECT_EQ(cm.capacity(), 0);
    cm.add(particles, 0, particles.size() - 1);
    EXPECT_EQ(cm.size(), particles.size());
    // EXPECT_EQ(cm.capacity(), cm.size()); // not true with default vector
}



TEST(CellMap, givesAccessToAllParticlesInACell)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, bucket_size, int, Point<int, 2>> cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    for (auto const& cell : patchbox)
    {
        auto& blist = cm.list_at(cell)->get();
        EXPECT_EQ(blist.size(), nppc);
        for (auto particleIndex : blist)
        {
            EXPECT_EQ(Point{particles[particleIndex].iCell}, cell);
        }
    }
}



TEST(CellMap, emptyLeavesCapacityButZeroSize)
{
    auto constexpr dim         = 2u;
    auto constexpr bucket_size = 100u;
    Box<int, 2> patchbox{{10, 20}, {25, 42}};
    CellMap<dim, bucket_size> cm;
    auto nppc      = 100u;
    auto particles = make_particles_in(patchbox, nppc);
    cm.add(particles);

    cm.empty();
    EXPECT_EQ(cm.size(), 0);
    EXPECT_TRUE(cm.is_empty());
    // EXPECT_EQ(cm.capacity(), patchbox.size() * nppc); // not true with default vector
}




// TEST(CellMap, trimMemory)
//{
//     auto constexpr dim         = 2u;
//     auto constexpr bucket_size = 100u;
//     Box<int, 2> patchbox{{10, 20}, {25, 42}};
//     CellMap<dim, bucket_size, int, Point<int, dim>> cm;
//     auto nppc      = 200u;
//     auto particles = make_particles_in(patchbox, nppc);
//     cm.add(particles);
//
//     EXPECT_EQ(cm.size(), patchbox.size() * nppc);
//     EXPECT_EQ(patchbox.size() * nppc, cm.capacity());
//
//     cm.empty();
//     EXPECT_EQ(cm.size(), 0);
//     EXPECT_EQ(cm.capacity(), patchbox.size() * nppc);
//     nppc      = 100;
//     particles = make_particles_in(patchbox, nppc);
//     cm.add(particles);
//     EXPECT_EQ(cm.size(), patchbox.size() * nppc);
//     EXPECT_EQ(patchbox.size() * nppc * 2, cm.capacity());
//     EXPECT_EQ(0.5, cm.used_mem_ratio());
//     cm.trim(1);
//     EXPECT_EQ(patchbox.size() * nppc * 2, cm.capacity());
//     cm.trim(0);
//     EXPECT_EQ(patchbox.size() * nppc * 1, cm.capacity());
// }


// TEST(CellMap, emptyTrimsWithMaxEmpty1IfCapacityExceeds3TimesSize)
//{
//     auto constexpr dim         = 2u;
//     auto constexpr bucket_size = 100u;
//     Box<int, 2> patchbox{{10, 20}, {25, 42}};
//     using cellmap_t = CellMap<dim, bucket_size, int, Point<int, 2>>;
//     cellmap_t cm;
//     auto nppc      = 600u;
//     auto particles = make_particles_in(patchbox, nppc);
//     cm.add(particles);
//     cm.empty();
//     nppc      = 100u;
//     particles = make_particles_in(patchbox, nppc);
//     cm.add(particles);
//     EXPECT_EQ(cm.capacity(), 600 * patchbox.size());
//     cm.empty();
//     EXPECT_EQ(cm.capacity(), particles.size() * 2);
// }



class CellMapExportFix : public ::testing::Test
{
public:
    CellMapExportFix()
    {
        particles = make_particles_in(patchbox, nppc);
        cm.add(particles);
    }

protected:
    static std::size_t constexpr dim         = 2;
    static std::size_t constexpr bucket_size = 100;
    static std::size_t constexpr nppc        = 100;
    using cellmap_t                          = CellMap<dim, bucket_size, int, Point<int, dim>>;
    cellmap_t cm;
    std::vector<Particle<dim>> particles;
    Box<int, dim> patchbox{{10, 20}, {25, 42}};
};


TEST_F(CellMapExportFix, exportItems)
{
    Box<int, 2> selectionBox{{14, 28}, {18, 37}};
    std::vector<Particle<dim>> selected;

    auto capa = cm.size(selectionBox);
    selected.reserve(capa);

    EXPECT_EQ(selected.size(), 0u);
    EXPECT_EQ(selected.capacity(), capa);

    cm.export_to(selectionBox, particles, selected);

    EXPECT_EQ(selected.size(), selected.capacity());
    for (auto const& p : selected)
    {
        EXPECT_TRUE(isIn(Point{p.iCell}, selectionBox));
    }
}

TEST_F(CellMapExportFix, exportWithTransform)
{
    Box<int, 2> selectionBox{{14, 28}, {18, 37}};
    std::vector<Particle<dim>> selected;
    auto capa = cm.size(selectionBox);
    selected.reserve(capa);
    EXPECT_EQ(selected.size(), 0u);
    EXPECT_EQ(selected.capacity(), capa);
    cm.export_to(selectionBox, particles, selected, [&](auto const& part) {
        auto copy{part};
        copy.iCell[0] += 100;
        return copy;
    });
    EXPECT_EQ(selected.size(), selected.capacity());
    auto offsetedSelectionBox{selectionBox};
    offsetedSelectionBox.lower[0] += 100;
    offsetedSelectionBox.upper[0] += 100;
    for (auto const& p : selected)
    {
        EXPECT_TRUE(isIn(Point{p.iCell}, offsetedSelectionBox));
    }
}


TEST_F(CellMapExportFix, exportWithPredicate)
{
    Box<int, 2> selectionBox{{14, 28}, {18, 37}};
    std::vector<Particle<dim>> selected;
    cm.export_to(particles, selected, [&](auto const& cell) { return isIn(cell, selectionBox); });
    for (auto const& p : selected)
    {
        EXPECT_TRUE(isIn(Point{p.iCell}, selectionBox));
    }
}


class CellMappedParticleBox : public ::testing::Test
{
public:
    CellMappedParticleBox()
        : patchBox{{10, 20, 30}, {25, 42, 54}}
        , ghostBox{patchBox}
        , outBox{patchBox}
    {
        ghostBox.grow(2);
        outBox.grow(4);
        particles = make_particles_in(outBox, nppc);
        cm.add(particles);
    }

    auto isInPatch()
    {
        return [this](auto const& cell) { return isIn(cell, patchBox); };
    }
    auto isInGhost()
    {
        return [this](auto const& cell) { return isIn(cell, ghostBox); };
    }

    auto isOut()
    {
        return [this](auto const& cell) { return isIn(cell, outBox) and !isIn(cell, ghostBox); };
    }

protected:
    static std::size_t constexpr dim         = 3;
    static std::size_t constexpr bucket_size = 100;
    static std::size_t constexpr nppc        = 1;
    using cellmap_t                          = CellMap<dim, bucket_size, int, Point<int, dim>>;
    Box<int, 3> patchBox;
    Box<int, 3> ghostBox;
    Box<int, 3> outBox;
    std::vector<Particle<dim>> particles;
    cellmap_t cm;
};




TEST_F(CellMappedParticleBox, trackParticle)
{
    EXPECT_EQ(cm.size(), particles.size());
    // pretends the particle change cell in x
    auto oldcell            = particles[200].iCell;
    particles[200].iCell[0] = 100;

    cm.update(particles, 200, oldcell);

    EXPECT_EQ(cm.size(), particles.size());

    auto blist = cm.list_at(particles[200].iCell);
    auto found = false;
    if (blist)
    {
        for (auto particleIndex : blist->get())
        {
            if (particleIndex == 200)
                found = true;
        }
    }
    EXPECT_TRUE(found);
}



TEST(CellMap, canReserveMemory)
{
    auto constexpr dim         = 3u;
    auto constexpr bucket_size = 100;
    using cellmap_t            = CellMap<dim, bucket_size, int, Point<int, dim>>;
    Box<int, 3> patchbox{{10, 20, 30}, {25, 42, 54}};
    cellmap_t cm;
    cm.reserve(patchbox);
    EXPECT_EQ(cm.nbr_cells(), patchbox.size());
}



TEST_F(CellMappedParticleBox, partitionsParticlesInPatchBox)
{
    EXPECT_EQ(cm.size(), particles.size());

    auto isInPatchBox = [&](auto const& cell) { return isIn(cell, patchBox); };
    auto allParts     = makeIndexRange(particles);
    auto inPatchRange = cm.partition(allParts, isInPatchBox);

    EXPECT_EQ(inPatchRange.size(), nppc * (patchBox.size()));
    EXPECT_EQ(particles.size() - inPatchRange.size(), nppc * (outBox.size() - patchBox.size()));

    for (std::size_t idx = inPatchRange.ibegin(); idx < inPatchRange.iend(); ++idx)
    {
        EXPECT_TRUE(isIn(Point{particles[idx].iCell}, patchBox));
    }
    for (std::size_t idx = inPatchRange.iend(); idx < particles.size(); ++idx)
    {
        EXPECT_FALSE(isIn(Point{particles[idx].iCell}, patchBox));
    }
}

TEST_F(CellMappedParticleBox, allButOneParticleSatisfyPredicate)
{
    EXPECT_EQ(cm.size(), particles.size());

    auto doNotTakeThatParticle = [&](auto const& cell) { return cell != patchBox.lower; };
    auto allParts              = makeIndexRange(particles);
    auto singleParticleRange   = cm.partition(allParts, doNotTakeThatParticle);

    EXPECT_EQ(singleParticleRange.size(), particles.size() - 1);
    EXPECT_EQ(singleParticleRange.ibegin(), 0);
    EXPECT_EQ(singleParticleRange.iend(), particles.size() - 1);

    for (std::size_t idx = singleParticleRange.ibegin(); idx < singleParticleRange.iend(); ++idx)
    {
        EXPECT_NE(Point{particles[idx].iCell}, patchBox.lower);
    }
}

TEST_F(CellMappedParticleBox, allButLastAlreadySatisfyPredicate)
{
    EXPECT_EQ(cm.size(), particles.size());

    auto lastParticleCell      = particles[particles.size() - 1].iCell;
    auto doNotTakeThatParticle = [&](auto const& cell) { return cell != lastParticleCell; };
    auto allParts              = makeIndexRange(particles);
    auto singleParticleRange   = cm.partition(allParts, doNotTakeThatParticle);

    EXPECT_EQ(singleParticleRange.size(), particles.size() - 1);
    EXPECT_EQ(singleParticleRange.ibegin(), 0);
    EXPECT_EQ(singleParticleRange.iend(), particles.size() - 1);

    for (std::size_t idx = singleParticleRange.ibegin(); idx < singleParticleRange.iend(); ++idx)
    {
        EXPECT_NE(Point{particles[idx].iCell}, lastParticleCell);
    }
}

TEST_F(CellMappedParticleBox, allOfTheParticlesSatisfyPredicate)
{
    EXPECT_EQ(cm.size(), particles.size());

    auto takeNoParticle    = [&](auto const& cell) { return true; };
    auto allParts          = makeIndexRange(particles);
    auto allParticlesRange = cm.partition(allParts, takeNoParticle);

    EXPECT_EQ(allParticlesRange.size(), particles.size());
    EXPECT_EQ(allParticlesRange.ibegin(), 0);
    EXPECT_EQ(allParticlesRange.iend(), particles.size());
}

TEST_F(CellMappedParticleBox, noneOfTheParticlesSatisfyPredicate)
{
    EXPECT_EQ(cm.size(), particles.size());

    auto takeNoParticle  = [&](auto const& cell) { return false; };
    auto allParts        = makeIndexRange(particles);
    auto noParticleRange = cm.partition(allParts, takeNoParticle);

    EXPECT_EQ(noParticleRange.size(), 0);
    EXPECT_EQ(noParticleRange.ibegin(), 0);
    EXPECT_EQ(noParticleRange.iend(), 0);
}

TEST_F(CellMappedParticleBox, rangeBasedPartition)
{
    auto partRange = makeRange(particles, 0, particles.size() / 2);
    auto inpatch   = cm.partition(partRange, isInPatch());

    // all particles in range before pivot should be in patchBox
    // but those in range after pivot should be outside

    for (std::size_t idx = inpatch.ibegin(); idx < inpatch.iend(); ++idx)
    {
        EXPECT_TRUE(isIn(Point{partRange.array()[idx].iCell}, patchBox));
    }
    for (std::size_t idx = inpatch.iend(); idx < partRange.iend(); ++idx)
    {
        EXPECT_FALSE(isIn(Point{partRange.array()[idx].iCell}, patchBox));
    }
}



TEST_F(CellMappedParticleBox, getPatchParticlesFromNonLeavingPartition)
{
    auto allParts = makeIndexRange(particles);

    // first get all particles still in ghost box
    // then from all those in ghostbox, take those still in patch
    auto inGhostBoxRange = cm.partition(allParts, isInGhost());
    auto inPatchRange    = cm.partition(inGhostBoxRange, isInPatch());

    EXPECT_EQ(inGhostBoxRange.size(), nppc * ghostBox.size());
    for (auto idx = inGhostBoxRange.ibegin(); idx < inGhostBoxRange.iend(); ++idx)
    {
        EXPECT_TRUE(isIn(Point{particles[idx].iCell}, ghostBox));
    }

    for (auto idx = inGhostBoxRange.iend(); idx < particles.size(); ++idx)
    {
        EXPECT_TRUE(isIn(Point{particles[idx].iCell}, outBox)
                    and !isIn(Point{particles[idx].iCell}, ghostBox));
    }

    for (auto idx = inPatchRange.ibegin(); idx < inPatchRange.iend(); ++idx)
    {
        EXPECT_TRUE(isIn(Point{particles[idx].iCell}, patchBox));
    }

    EXPECT_EQ(0, inPatchRange.ibegin());
    EXPECT_EQ(nppc * patchBox.size(), inPatchRange.size());
}


TEST_F(CellMappedParticleBox, eraseOutOfPatchRange)
{
    auto allParts = makeIndexRange(particles);
    auto inpatch  = cm.partition(allParts, isInPatch());
    auto toErase  = makeRange(particles, inpatch.iend(), particles.size());
    cm.erase(toErase);

    for (auto const& part : particles)
        EXPECT_TRUE(isIn(Point{part.iCell}, patchBox));

    EXPECT_EQ(particles.size(), patchBox.size() * nppc);
}




#if 0
// keep it warm here maybe useful in the future if/when sorting is needed
TEST(CellMap, sortArray)
{
    auto constexpr dim         = 3u;
    auto constexpr bucket_size = 100;
    using cellmap_t            = CellMap<dim, bucket_size, int, Point<int, dim>>;
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
#endif

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

