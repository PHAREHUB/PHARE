// The ion updater on every patch of a level, with each patch's selection boxes made as the
// solver makes them. Particles leaving a patch must either stay in its domain or become patch
// ghosts, for the messengers to hand them to the patch that owns their cell: none may be lost.

#include "phare_core.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep
#include "simulator/simulator_def.hpp"

#include "core/numerics/ion_updater/ion_updater.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "tests/amr/test_level_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include <memory>
#include <vector>


using namespace PHARE;

template<std::size_t dim_, std::size_t interp_>
struct UpdaterOnLevel
{
    auto static constexpr dim    = dim_;
    auto static constexpr interp = interp_;
    SimOpts static constexpr opts{dim, interp};

    using PHARE_Types     = core::PHARE_Types<opts>;
    using GridLayout_t    = PHARE_Types::Hybrid::GridLayout_t;
    using ParticleArray_t = PHARE_Types::Hybrid::ParticleArray_t;
    using Particle_t      = ParticleArray_t::value_type;
    using Box_t           = core::Box<int, dim>;
    using Electromag_t    = core::UsableElectromag<dim>;
    using Ions_t          = core::UsableIons_t<ParticleArray_t, interp>;
    using IonUpdater_t    = core::IonUpdater<Ions_t, Electromag_t, GridLayout_t>;
    using Boxing_t        = core::UpdaterSelectionBoxing<IonUpdater_t, GridLayout_t>;

    // as SolverPPC::make_boxes, E = B = 0 so particles move in straight lines
    struct Patch
    {
        Patch(SAMRAI::hier::Patch const& patch, SAMRAI::hier::HierarchyNeighbors const& neighbors)
            : layout{amr::layoutFromPatch<GridLayout_t>(patch)}
            , boxing{layout, amr::makeNonLevelGhostBoxFor<GridLayout_t>(patch, neighbors),
                     amr::makeNonOwnedBoxesFor<dim>(patch, neighbors)}
            , em{layout}
            , ions{layout, "protons"}
        {
            em.E.zero();
            em.B.zero();
        }

        auto& domain() { return ions.populations[0].particles.domain_particles; }
        auto& patchGhosts() { return ions.populations[0].particles.patch_ghost_particles; }

        GridLayout_t layout;
        Boxing_t boxing;
        Electromag_t em;
        Ions_t ions; // not movable: buffers point into it
    };

    // one particle starting in `cell` of the patch of box `from`, moving by `cells` in one push
    struct Move
    {
        Box_t from;
        std::array<int, dim> cell;
        std::array<double, dim> cells;
    };

    UpdaterOnLevel(Box_t const& domain, std::vector<Box_t> const& boxes)
        : boxes{boxes}
        , hierarchy{test::hierarchyWithLevel0(domain, boxes)}
    {
        auto const neighbors = amr::makeSameLevelNeighbors<GridLayout_t>(*hierarchy, 0);
        for (auto const& patch : *hierarchy->getPatchLevel(0))
            patches.emplace_back(std::make_unique<Patch>(*patch, neighbors));
    }

    Patch& patchOf(Box_t const& box)
    {
        for (auto& patch : patches)
            if (patch->layout.AMRBox() == box)
                return *patch;
        throw std::runtime_error("no patch for box");
    }

    // particle velocities in cells per time step, the mesh size is 1
    void add(Move const& move)
    {
        std::array<double, dim> delta;
        delta.fill(.5);
        std::array<double, 3> v{0, 0, 0};
        for (std::size_t i = 0; i < dim; ++i)
            v[i] = move.cells[i] / dt;
        patchOf(move.from).domain().push_back(Particle_t{/*weight*/ 1., move.cell, delta, v});
    }

    void update()
    {
        initializer::PHAREDict dict;
        dict["pusher"]["name"] = std::string{"modified_boris"};
        IonUpdater_t ionUpdater{dict};
        for (auto& patch : patches)
            ionUpdater.updatePopulations(patch->ions, patch->em, patch->boxing, dt,
                                         core::UpdaterMode::all);
    }

    // every particle is in the domain of its patch or a patch ghost inside another patch
    void checkNoneLost(std::size_t const nbrParticles)
    {
        std::size_t found = 0;
        for (auto& patch : patches)
        {
            auto const box = patch->layout.AMRBox();
            for (auto const& particle : patch->domain())
                EXPECT_TRUE(isIn(core::Point{particle.iCell}, box));
            for (auto const& particle : patch->patchGhosts())
                EXPECT_TRUE(std::any_of(boxes.begin(), boxes.end(),
                                        [&](auto const& other) {
                                            return !(other == box)
                                                   and isIn(core::Point{particle.iCell}, other);
                                        }))
                    << "patch ghost " << core::Point{particle.iCell} << " of " << box
                    << " is in no other patch";
            found += patch->domain().size() + patch->patchGhosts().size();
        }
        EXPECT_EQ(found, nbrParticles) << nbrParticles - found << " particle(s) lost";
    }

    double const dt = 1;
    std::vector<Box_t> const boxes;
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy;
    std::vector<std::unique_ptr<Patch>> patches;
};


// C touches A only at A's upper right corner, D (right of A) is one cell taller than A:
//   B [0..9]x[10..20]   C [10..19]x[11..20]
//   A [0..9]x[0..9]     D [10..19]x[0..10]
template<std::size_t interp>
void checkCornerNeighbourParticlesAreNotLost()
{
    using Level = UpdaterOnLevel<2, interp>;
    using Box_t = Level::Box_t;
    Box_t const A{{0, 0}, {9, 9}}, B{{0, 10}, {9, 20}}, C{{10, 11}, {19, 20}}, D{{10, 0}, {19, 10}};

    Level level{{{0, 0}, {19, 20}}, {A, B, D, C}};
    std::vector<typename Level::Move> const moves{
        {A, {9, 9}, {1, 2}},     // into C (10,11), only at A's corner
        {C, {10, 11}, {-1, -2}}, // into A (9,9), only at C's corner
        {A, {9, 9}, {1, 0}},     // into D (10,9), adjacent
        {C, {10, 15}, {-1, 0}},  // into B (9,15), adjacent
    };
    for (auto const& move : moves)
        level.add(move);

    level.update();
    level.checkNoneLost(moves.size());
}

// moving 2 cells in a push needs a particle ghost width of 2: interp 2 and 3
TEST(IonUpdaterOnLevel, cornerNeighbourParticlesAreNotLostInterp2)
{
    checkCornerNeighbourParticlesAreNotLost<2>();
}
TEST(IonUpdaterOnLevel, cornerNeighbourParticlesAreNotLostInterp3)
{
    checkCornerNeighbourParticlesAreNotLost<3>();
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    int testResult = RUN_ALL_TESTS();

    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
