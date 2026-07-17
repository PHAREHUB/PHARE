
#include "phare_core.hpp"

#include "core/utilities/types.hpp"
#include "core/numerics/ion_updater/ion_updaters.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/numerics/interpolator/interpolating.hpp"
#include "core/numerics/interpolator/interpolating.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_converter.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"
#include "tests/amr/resources_manager/test_resources_manager_fixtures.hpp"

#include "gtest/gtest.h"

#include <unordered_map>


namespace PHARE::core
{

using Param  = std::vector<double> const&;
using Return = std::shared_ptr<PHARE::core::Span<double>>;

Return density(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}

Return vx(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return vy(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return vz(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return vthx(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), .1);
}

Return vthy(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), .1);
}

Return vthz(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), .1);
}

Return bx(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return by(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}

Return bz(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 0);
}



std::size_t nbrPartPerCell = 1000;

using InitFunctionT = PHARE::initializer::InitFunction<1>;

PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;

    dict["simulation"]["algo"]["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};

    dict["ions"]["nbrPopulations"]                          = std::size_t{2};
    dict["ions"]["pop0"]["name"]                            = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]                            = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["name"]    = std::string{"maxwellian"};
    dict["ions"]["pop0"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(density);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(vx);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(vy);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(vz);


    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(vthx);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(vthy);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(vthz);


    dict["ions"]["pop0"]["particle_initializer"]["nbr_part_per_cell"]
        = static_cast<int>(nbrPartPerCell);
    dict["ions"]["pop0"]["particle_initializer"]["charge"] = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["basis"]  = std::string{"cartesian"};

    dict["ions"]["pop1"]["name"]                            = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"]                            = 1.;
    dict["ions"]["pop1"]["particle_initializer"]["name"]    = std::string{"maxwellian"};
    dict["ions"]["pop1"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(density);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(vx);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(vy);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(vz);


    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(vthx);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(vthy);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(vthz);


    dict["ions"]["pop1"]["particle_initializer"]["nbr_part_per_cell"]
        = static_cast<int>(nbrPartPerCell);
    dict["ions"]["pop1"]["particle_initializer"]["charge"] = 1.;
    dict["ions"]["pop1"]["particle_initializer"]["basis"]  = std::string{"cartesian"};

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<InitFunctionT>(bx);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<InitFunctionT>(by);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<InitFunctionT>(bz);

    return dict;
}
static auto init_dict = createDict();




template<auto opts_>
struct TestParam
{
    auto static constexpr opts         = opts_;
    static constexpr auto dimension    = opts.dimension;
    static constexpr auto interp_order = opts.interp_order;
};

template<auto layout_mode, auto opts>
struct updater_test_bits;

template<auto opts>
struct updater_test_bits<core::LayoutMode::AoSTS, opts>
{
    auto static constexpr atomic_ops = false;
    using PHARETypes                 = PHARE::core::PHARE_Types<opts>;
    using Ions                       = PHARETypes::Ions_t;
    using Electromag                 = PHARETypes::Electromag_t;
    using GridLayout                 = PHARETypes::GridLayout_t;
    using ParticleArray              = PHARETypes::ParticleArray_t;
    using IonUpdater = core::IonUpdaterImplResolver<Ions, Electromag, GridLayout>::IonUpdater_t;
    using Boxing_t   = PHARE::core::UpdaterSelectionBoxing<GridLayout>;
};
template<auto opts>
struct updater_test_bits<core::LayoutMode::AoSMapped, opts>
{
    auto static constexpr atomic_ops = false;
    using PHARETypes                 = PHARE::core::PHARE_Types<opts>;
    using Ions                       = PHARETypes::Ions_t;
    using Electromag                 = PHARETypes::Electromag_t;
    using GridLayout                 = PHARETypes::GridLayout_t;
    using ParticleArray              = PHARETypes::ParticleArray_t;
    using IonUpdater = core::IonUpdaterImplResolver<Ions, Electromag, GridLayout>::IonUpdater_t;
    using Selector_t = IonUpdater::Pusher::ParticleSelector;
    using Boxing_t   = PHARE::core::UpdaterCellMapSelectionBoxing<Selector_t, GridLayout>;
};

template<auto layout_mode, auto opts>
struct updater_test_bits
{
    // finish if needed
};

template<typename TestParam_t>
struct IonUpdaterTest : public ::testing::Test,
                        public updater_test_bits<TestParam_t::opts.layout_mode, TestParam_t::opts>
{
    auto static constexpr opts         = TestParam_t::opts;
    static constexpr auto dim          = opts.dimension;
    static constexpr auto interp_order = opts.interp_order;
    using PHARETypes                   = PHARE::core::PHARE_Types<opts>;
    using Ions                         = PHARETypes::Ions_t;
    using Electromag                   = PHARETypes::Electromag_t;
    using GridLayout                   = PHARETypes::GridLayout_t;
    using ParticleArray                = PHARETypes::ParticleArray_t;
    using ParticleInitializerFactory   = PHARETypes::ParticleInitializerFactory_t;

    using basics       = updater_test_bits<TestParam_t::opts.layout_mode, TestParam_t::opts>;
    using IonUpdater   = basics::IonUpdater;
    using Boxing_t     = basics::Boxing_t;
    using Interpolator = Interpolating<ParticleArray, interp_order, basics::atomic_ops>;

    using UsableElectromag_t = UsableElectromag<GridLayout, opts.alloc_mode, opts.layout_mode>;

    Interpolator interpolator;

    double const dt{0.01};

    // grid configuration
    GridLayout const layout{{0.1}, {100u}, {{0.}}};


    Boxing_t const boxing{layout, {layout.AMRBox()}};
    std::unordered_map<std::string, Boxing_t> const levelBoxing{{"patch_id", boxing}};

    static_assert(std::is_same_v<Electromag, typename UsableElectromag_t::Super>);
    UsableElectromag_t EM{layout, init_dict["electromag"]};

    UsableIons_t<ParticleArray, interp_order> ions{layout, init_dict["ions"]};

    IonUpdaterTest()
    {
        // ok all resources pointers are set to buffers
        // now let's initialize Electromag fields to user input functions
        // and ion population particles to user supplied moments

        assert(no_nans(ions.velocity()(Component::X)));
        EM.initialize(layout);
        for (auto& pop : ions)
        {
            auto const& info         = pop.particleInitializerInfo();
            auto particleInitializer = ParticleInitializerFactory::create(info);
            particleInitializer->loadParticles(pop.domainParticles(), layout);
            EXPECT_GT(pop.domainParticles().size(), 0ull);
            // particle_array_domain_is_valid(pop.domainParticles(), layout.AMRBox());
        }
        assert(no_nans(ions.velocity()(Component::X)));

        // now all domain particles are loaded we need to manually insert
        // ghost particles (this is in reality SAMRAI's job)
        // these are needed if we want all used nodes to be complete


        // in 1D we assume left border is touching the level border
        // and right is touching another patch
        // so on the left no patchGhost but levelGhost(and old and new)
        // on the right no levelGhost but patchGhosts


        for (auto& pop : ions)
        {
            if constexpr (dim == 1)
            {
                int firstPhysCell = layout.physicalStartIndex(QtyCentering::dual, Direction::X);
                int lastPhysCell  = layout.physicalEndIndex(QtyCentering::dual, Direction::X);
                auto firstAMRCell = layout.localToAMR(Point{firstPhysCell});
                auto lastAMRCell  = layout.localToAMR(Point{lastPhysCell});

                // we need to put levelGhost particles in the cell just to the
                // left of the first cell. In reality these particles should
                // come from splitting particles of the next coarser level
                // in this test we just copy those of the first cell
                // we also assume levelGhostOld and New are the same particles
                // for simplicity

                auto& domainPart        = pop.domainParticles();
                auto& levelGhostPartOld = pop.levelGhostParticlesOld();
                auto& levelGhostPartNew = pop.levelGhostParticlesNew();
                auto& levelGhostPart    = pop.levelGhostParticles();
                auto& patchGhostPart    = pop.patchGhostParticles();



                // copies need to be put in the ghost cell
                // we have copied particles be now their iCell needs to be udpated
                // our choice is :
                //
                // first order:
                //
                //   ghost| domain...
                // [-----]|[-----][-----][-----][-----][-----]
                //     ^      v
                //     |      |
                //     -------|
                //
                // second and third order:

                //   ghost        | domain...
                // [-----]|[-----][-----][-----][-----][-----][-----]
                //     ^      ^       v     v
                //     |      |       |     |
                //     -------|-------|     |
                //            ---------------

                auto const per_particles = [&](auto& parts) {
                    for (auto const& part : parts)
                    {
                        if constexpr (interp_order == 2 or interp_order == 3)
                        {
                            if (part.iCell()[0] == firstAMRCell[0]
                                or part.iCell()[0] == firstAMRCell[0] + 1)
                            {
                                auto p{part};
                                p.iCell()[0] -= 2;
                                levelGhostPartOld.push_back(p);
                            }
                        }
                        else if constexpr (interp_order == 1)
                        {
                            if (part.iCell()[0] == firstAMRCell[0])
                            {
                                auto p{part};
                                p.iCell()[0] -= 1;
                                levelGhostPartOld.push_back(p);
                            }
                        }
                    }
                };

                if constexpr (opts.layout_mode == LayoutMode::AoSTS)
                {
                    for (auto& tile : domainPart())
                        per_particles(tile());
                }
                else
                {
                    per_particles(domainPart);
                }

                ParticleArrayService::template sync<0, ParticleType::Ghost>(levelGhostPartOld);
                levelGhostPartNew = levelGhostPartOld;
                levelGhostPart    = levelGhostPartOld;

                EXPECT_EQ(pop.domainParticles().size(), 100 * nbrPartPerCell);
                EXPECT_EQ(levelGhostPartOld.size(),
                          GridLayout::nbrParticleGhosts() * nbrPartPerCell);
                EXPECT_EQ(patchGhostPart.size(), 0);


            } // end 1D
        } // end pop loop

        assert(no_nans(ions.velocity()(Component::X)));
        PHARE::core::depositParticles(ions, layout, interpolator, PHARE::core::DomainDeposit{});
        PHARE::core::depositParticles(ions, layout, interpolator, PHARE::core::LevelGhostDeposit{});
        assert(no_nans(ions.velocity()(Component::X)));

        ions.computeChargeDensity();
        assert(no_nans(ions.velocity()(Component::X)));
        ions.computeBulkVelocity();
        assert(no_nans(ions.velocity()(Component::X)));
    } // end Ctor



    using Patch = test::HybridPatch<Ions, Electromag>;

    Patch as_patch() { return Patch{layout, *ions, *EM}; }

    auto update(auto& ionUpdater, auto const mode)
    {
        assert(no_nans(ions.velocity()(Component::X)));
        if constexpr (ParticleArray::layout_mode == AoSMapped)
            ionUpdater.updatePopulations(*ions, *EM, boxing, dt, mode);

        else
        {
            std::vector<Patch> patches{as_patch()};

            auto const quantities = [&](int i) {
                return std::forward_as_tuple(patches[0].model.state.ions,
                                             patches[i].model.state.electromag);
            };

            assert(ions.massDensity().data() == patches[0].model.state.ions.massDensity().data());
            auto accessor = test::make_model_level_accessor(patches, quantities);
            ionUpdater.updatePopulations(accessor, levelBoxing, dt, mode);
        }
        assert(no_nans(ions.velocity()(Component::X)));
    }

    void fillIonsMomentsGhosts()
    {
        for (auto& pop : this->ions)
        {
            double alpha = 0.5;
            interpolator(pop.levelGhostParticlesNew(), pop.particleDensity(), pop.chargeDensity(),
                         pop.flux(), layout,
                         /*coef = */ alpha);

            interpolator(pop.levelGhostParticlesOld(), pop.particleDensity(), pop.chargeDensity(),
                         pop.flux(), layout,
                         /*coef = */ (1. - alpha));
        }
    }



    void checkMomentsHaveEvolved(auto const& ionsBufferCpy)
    {
        auto& populations = this->ions.getRunTimeResourcesViewList();

        auto const& protonParticleDensity = reduce(populations[0].particleDensity());
        auto const& protonChargeDensity   = reduce(populations[0].chargeDensity());
        auto const& protonFx              = reduce(populations[0].flux()(Component::X));
        auto const& protonFy              = reduce(populations[0].flux()(Component::Y));
        auto const& protonFz              = reduce(populations[0].flux()(Component::Z));

        auto const& alphaParticleDensity = reduce(populations[1].particleDensity());
        auto const& alphaChargeDensity   = reduce(populations[1].chargeDensity());
        auto const& alphaFx              = reduce(populations[1].flux()(Component::X));
        auto const& alphaFy              = reduce(populations[1].flux()(Component::Y));
        auto const& alphaFz              = reduce(populations[1].flux()(Component::Z));


        auto ix0 = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
        auto ix1 = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

        auto nonZero = [&](auto const& field, std::string const type) {
            auto sum = 0.;
            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                sum += std::abs(field(ix));
            }
            EXPECT_GT(sum, 0.);
            if (sum == 0)
                std::cout << "nonZero failed for " << type << " : " << field.name() << "\n";
        };

        auto check = [&](auto const& newField, auto const& og) {
            // assert(no_nans(newField));
            assert(no_nans(og));
            auto const& originalField = reduce(og);
            nonZero(newField, newField.name() + ":new");
            // nonZero(originalField, "originalField");
            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                auto evolution = std::abs(newField(ix) - originalField(ix));
                //  should check that moments are still compatible with user inputs also
                EXPECT_GT(evolution, 0.0);
                if (evolution <= 0.0)
                    std::cout << newField.name() + ": after update : " << newField(ix)
                              << " before update : " << originalField(ix)
                              << " evolution : " << evolution << " ix : " << ix << "\n";
            }
        };


        check(protonParticleDensity, ionsBufferCpy[0].particleDensity());
        check(protonChargeDensity, ionsBufferCpy[0].chargeDensity());
        check(protonFx, ionsBufferCpy[0].flux()(Component::X));
        check(protonFy, ionsBufferCpy[0].flux()(Component::Y));
        check(protonFz, ionsBufferCpy[0].flux()(Component::Z));

        check(alphaParticleDensity, ionsBufferCpy[1].particleDensity());
        check(alphaChargeDensity, ionsBufferCpy[1].chargeDensity());
        check(alphaFx, ionsBufferCpy[1].flux()(Component::X));
        check(alphaFy, ionsBufferCpy[1].flux()(Component::Y));
        check(alphaFz, ionsBufferCpy[1].flux()(Component::Z));

        assert(no_nans(ionsBufferCpy.velocity()(Component::X)));
        check(reduce(ions.velocity()(Component::X)), ionsBufferCpy.velocity()(Component::X));
        check(reduce(ions.velocity()(Component::Y)), ionsBufferCpy.velocity()(Component::Y));
        check(reduce(ions.velocity()(Component::Z)), ionsBufferCpy.velocity()(Component::Z));
    }



    void checkDensityIsAsPrescribed()
    {
        auto ix0 = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
        auto ix1 = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

        auto check = [&](auto const& density, auto const& function) {
            std::vector<std::size_t> ixes;
            std::vector<double> x;

            // We do not use the primal box as the last primal is incomplete.
            // This is because level ghosts are only interpolated if they enter the domain
            for (auto const [amr_idx, lcl_idx] : this->layout.amr_lcl_idx())
            {
                ixes.emplace_back(lcl_idx[0]);
                x.emplace_back(layout.cellCenteredCoordinates(amr_idx)[0]);
            }

            auto functionXPtr = function(x); // keep alive
            EXPECT_EQ(functionXPtr->size(), (ix1 - ix0));

            auto& functionX = *functionXPtr;

            for (std::size_t i = 0; i < functionX.size(); i++)
            {
                auto ix   = ixes[i];
                auto diff = std::abs(density(ix) - functionX[i]);

                EXPECT_GE(0.07, diff);

                if (diff >= 0.07)
                    std::cout << i << " actual : " << density(ix)
                              << " prescribed : " << functionX[i] << " diff : " << diff
                              << " ix : " << ix << " max: " << functionX.size() << "\n";
            }
        };

        auto& populations           = this->ions.getRunTimeResourcesViewList();
        auto& protonParticleDensity = populations[0].particleDensity();
        auto& alphaParticleDensity  = populations[1].particleDensity();

        check(reduce(protonParticleDensity), density);
        check(reduce(alphaParticleDensity), density);
    }
};


// clang-format off
using Permutations = ::testing::Types<
    TestParam<SimOpts{1, 1}>
  , TestParam<SimOpts{1, 2}>
  , TestParam<SimOpts{1, 3}>
PHARE_WITH_MKN_GPU(
  , TestParam<SimOpts{.dimension = 1, .interp_order = 1, .layout_mode = LayoutMode::AoSTS}>
)
>;
// clang-format on
TYPED_TEST_SUITE(IonUpdaterTest, Permutations, );




TYPED_TEST(IonUpdaterTest, ionUpdaterTakesPusherParamsFromPHAREDictAtConstruction)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{};
}


// the following 3 tests are testing the fixture is well configured.


TYPED_TEST(IonUpdaterTest, loadsDomainPatchAndLevelGhostParticles)
{
    auto check = [this](std::size_t nbrGhostCells, auto& pop) {
        EXPECT_EQ(this->layout.nbrCells()[0] * nbrPartPerCell, pop.domainParticles().size());
        EXPECT_EQ(0, pop.patchGhostParticles().size());
        EXPECT_EQ(nbrGhostCells * nbrPartPerCell, pop.levelGhostParticlesOld().size());
        EXPECT_EQ(nbrGhostCells * nbrPartPerCell, pop.levelGhostParticlesNew().size());
        EXPECT_EQ(nbrGhostCells * nbrPartPerCell, pop.levelGhostParticles().size());
    };


    if constexpr (TypeParam::dimension == 1)
    {
        for (auto& pop : this->ions)
        {
            if constexpr (TypeParam::interp_order == 1)
            {
                check(1, pop);
            }
            else if constexpr (TypeParam::interp_order == 2 or TypeParam::interp_order == 3)
            {
                check(2, pop);
            }
        }
    }
}




TYPED_TEST(IonUpdaterTest, loadsLevelGhostParticlesOnLeftGhostArea)
{
    int firstPhysCell = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto firstAMRCell = this->layout.localToAMR(Point{firstPhysCell});

    if constexpr (TypeParam::dimension == 1)
    {
        for (auto& pop : this->ions)
        {
            auto const ghost = convert_to<LayoutMode::AoS>(pop.levelGhostParticles(), this->layout);

            if constexpr (TypeParam::interp_order == 1)
            {
                for (auto it = ghost.begin(); it != ghost.end(); ++it)
                {
                    EXPECT_EQ(firstAMRCell[0] - 1, it.iCell()[0]);
                }
            }
            else if constexpr (TypeParam::interp_order == 2 or TypeParam::interp_order == 3)
            {
                auto copy                 = ghost;
                auto firstInOuterMostCell = std::partition(
                    std::begin(copy), std::end(copy), [&firstAMRCell](auto const& particle) {
                        return particle.iCell()[0] == firstAMRCell[0] - 1;
                    });
                EXPECT_EQ(nbrPartPerCell, std::distance(std::begin(copy), firstInOuterMostCell));
                EXPECT_EQ(nbrPartPerCell, std::distance(firstInOuterMostCell, std::end(copy)));
            }
        }
    }
}




// start of PHARE TESTS



TYPED_TEST(IonUpdaterTest, particlesUntouchedInMomentOnlyMode)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{};

    auto ionsBufferCpy = this->ions;

    this->update(ionUpdater, UpdaterMode::domain_only);

    this->fillIonsMomentsGhosts();

    ionUpdater.updateIons(this->ions);

    auto& populations = this->ions.getRunTimeResourcesViewList();

    auto checkIsUnTouched = [&](auto const& og, auto const& cp) {
        auto const original = convert_to<LayoutMode::AoS>(og, this->layout);
        auto const cpy      = convert_to<LayoutMode::AoS>(cp, this->layout);

        // no particles should have moved, so none should have left the domain
        EXPECT_EQ(cpy.size(), original.size());
        for (std::size_t iPart = 0; iPart < original.size(); ++iPart)
        {
            EXPECT_EQ(cpy.iCell(iPart)[0], original.iCell(iPart)[0]);

            EXPECT_DOUBLE_EQ(cpy.delta(iPart)[0], original.delta(iPart)[0]);
            for (std::size_t iDir = 0; iDir < 3; ++iDir)
            {
                EXPECT_DOUBLE_EQ(cpy.v(iPart)[iDir], original.v(iPart)[iDir]);
            }
        }
    };

    auto& cpy_pops = ionsBufferCpy.getRunTimeResourcesViewList();

    checkIsUnTouched(populations[0].levelGhostParticles(), cpy_pops[0].levelGhostParticles());
    checkIsUnTouched(populations[0].levelGhostParticlesOld(), cpy_pops[0].levelGhostParticlesOld());
    checkIsUnTouched(populations[0].levelGhostParticlesNew(), cpy_pops[0].levelGhostParticlesNew());

    checkIsUnTouched(populations[1].levelGhostParticles(), cpy_pops[1].levelGhostParticles());
    checkIsUnTouched(populations[1].levelGhostParticlesOld(), cpy_pops[1].levelGhostParticlesOld());
    checkIsUnTouched(populations[1].levelGhostParticlesNew(), cpy_pops[1].levelGhostParticlesNew());
}




// TYPED_TEST(IonUpdaterTest, particlesAreChangedInParticlesAndMomentsMode)
//{
//    typename IonUpdaterTest<TypeParam>::IonUpdater
//    ionUpdater{init_dict["simulation"]["pusher"]};
//
//    IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};
//
//    ionUpdater.updatePopulations(this->ions, this->EM, this->boxing, this->dt,
//                                 UpdaterMode::particles_and_moments);
//
//    this->fillIonsMomentsGhosts();
//
//    ionUpdater.updateIons(this->ions);
//
//    auto& populations = this->ions.getRunTimeResourcesViewList();
//
//    EXPECT_NE(ionsBufferCpy.protonDomain.size(), populations[0].domainParticles().size());
//    EXPECT_NE(ionsBufferCpy.alphaDomain.size(), populations[1].domainParticles().size());
//
//    // cannot think of anything else to check than checking that the number of particles
//    // in the domain have changed after them having been pushed.
//}



TYPED_TEST(IonUpdaterTest, momentsAreChangedInParticlesAndMomentsMode)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{};

    auto ionsBufferCpy = this->ions;

    assert(ionsBufferCpy.massDensity().data() != this->ions.massDensity().data());

    assert(no_nans(this->ions.velocity()(Component::X)));
    this->update(ionUpdater, UpdaterMode::all);
    assert(no_nans(this->ions.velocity()(Component::X)));

    this->fillIonsMomentsGhosts();
    assert(no_nans(this->ions.velocity()(Component::X)));

    ionUpdater.updateIons(this->ions);

    // assert(no_nans(this->ions.velocity()(Component::X)));

    this->checkMomentsHaveEvolved(ionsBufferCpy);
    this->checkDensityIsAsPrescribed();
}




TYPED_TEST(IonUpdaterTest, momentsAreChangedInMomentsOnlyMode)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{};

    auto ionsBufferCpy = this->ions;

    assert(ionsBufferCpy.massDensity().data() != this->ions.massDensity().data());

    this->update(ionUpdater, UpdaterMode::domain_only);

    this->fillIonsMomentsGhosts();

    ionUpdater.updateIons(this->ions);

    this->checkMomentsHaveEvolved(ionsBufferCpy);
    this->checkDensityIsAsPrescribed();
}



TYPED_TEST(IonUpdaterTest, thatNoNaNsExistOnPhysicalNodesMoments)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater ionUpdater{};

    this->update(ionUpdater, UpdaterMode::domain_only);

    this->fillIonsMomentsGhosts();

    ionUpdater.updateIons(this->ions);

    auto ix0 = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto ix1 = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    for (auto& pop : this->ions)
    {
        for (auto ix = ix0; ix <= ix1; ++ix)
        {
            auto const& density = reduce(pop.particleDensity());
            auto const& flux    = pop.flux();

            auto const& fx = reduce(flux(Component::X));
            auto const& fy = reduce(flux(Component::Y));
            auto const& fz = reduce(flux(Component::Z));

            EXPECT_FALSE(std::isnan(density(ix)));
            EXPECT_FALSE(std::isnan(fx(ix)));
            EXPECT_FALSE(std::isnan(fy(ix)));
            EXPECT_FALSE(std::isnan(fz(ix)));
        }
    }
}


template<std::size_t dim>
void pseudo_randomize_deltas(auto& parts)
{
    std::size_t k = 0;

    per_particle(parts, [&](Particle<dim>& p) {
        for (std::uint16_t i = 0; i < dim; ++i)
        {
            double dir = k % 2 == 0 ? -1 : 1;
            p.delta()[i] += dir * .03 * k;
            assert(p.delta()[i] > 0 and p.delta()[i] < 1);
            if (k == 9)
                k = 0;
            ++k;
        }
    });
}

template<typename GridLayout>
void populate_particles(auto& ions, GridLayout& layout)
{
    auto constexpr dim = GridLayout::dimension;
    Particle<dim> const particle{1.0 / nbrPartPerCell, 1, ConstArray<int, dim>(),
                                 ConstArray<double, dim>(.5), ConstArray<double, 3>(1)};
    auto const ghost_box = grow(layout.AMRBox(), 1);

    auto const shift_particle = [](auto p, auto const bix) {
        p.iCell() = bix;
        return p;
    };

    for (auto& pop : ions)
    {
        pop.domainParticles().clear();
        pop.levelGhostParticles().clear();
        pop.patchGhostParticles().clear();

        pop.domainParticles().check();
        pop.levelGhostParticles().check();

        for (auto const& bix : layout.AMRBox())
            for (std::size_t i = 0; i < nbrPartPerCell; ++i)
                pop.domainParticles().emplace_back(shift_particle(particle, *bix));

        for (auto const& bix : ghost_box)
            if (not isIn(bix, layout.AMRBox()))
                for (std::size_t i = 0; i < nbrPartPerCell; ++i)
                    pop.levelGhostParticles().emplace_back(shift_particle(particle, *bix));

        per_particle(pop.levelGhostParticles(),
                     [&](Particle<dim>& p) { EXPECT_TRUE(not isIn(p, layout.AMRBox())); });

        pseudo_randomize_deltas<dim>(pop.domainParticles());
        pseudo_randomize_deltas<dim>(pop.levelGhostParticles());

        pop.domainParticles().check();
        pop.levelGhostParticles().check();
    }
}




// TYPED_TEST(IonUpdaterTest, hardCodedRegressionTest0)
// {
//     using IonUpdater      = TestFixture::IonUpdater;
//     auto constexpr dim    = TestFixture::dim;
//     auto constexpr interp = TestFixture::interp_order;

//     std::array<double, 3> expectations
//         = {199.99999999999889, 199.99999999999659, 199.99999999999673};

//     populate_particles(this->ions, this->layout);

//     for (auto& pop : this->ions)
//     {
//         EXPECT_EQ(pop.domainParticles().size(), 100 * nbrPartPerCell);
//         EXPECT_EQ(pop.levelGhostParticles().size(), 2 * nbrPartPerCell);
//         EXPECT_EQ(pop.patchGhostParticles().size(), 0);
//     }

//     IonUpdater ionUpdater{};
//     this->update(ionUpdater, UpdaterMode::domain_only);

//     std::int64_t icellSum = 0;
//     double densitySum = 0, fluxXSum = 0, fluxYSum = 0, fluxZSum = 0;
//     for (auto& pop : this->ions)
//     {
//         densitySum += PHARE::core::sum(reduce(pop.density()));
//         fluxXSum += PHARE::core::sum(reduce(pop.flux()[0]));
//         fluxYSum += PHARE::core::sum(reduce(pop.flux()[1]));
//         fluxZSum += PHARE::core::sum(reduce(pop.flux()[2]));
//         per_particle(pop.domainParticles(),
//                      [&](Particle<dim>& p) { icellSum += PHARE::core::sum(p.iCell()); });
//     }

//     EXPECT_EQ(icellSum, 9900000);
//     EXPECT_DOUBLE_EQ(densitySum, expectations[interp - 1]);
//     EXPECT_DOUBLE_EQ(fluxXSum, expectations[interp - 1]);
//     EXPECT_DOUBLE_EQ(fluxYSum, expectations[interp - 1]);
//     EXPECT_DOUBLE_EQ(fluxZSum, expectations[interp - 1]);
// }



// TYPED_TEST(IonUpdaterTest, hardCodedRegressionTest1)
// {
//     using IonUpdater      = TestFixture::IonUpdater;
//     auto constexpr dim    = TestFixture::dim;
//     auto constexpr interp = TestFixture::interp_order;

//     std::array<double, 3> expectations
//         = {199.99999999999889, 199.99999999999659, 199.99999999999673};

//     populate_particles(this->ions, this->layout);

//     double deltaSum = 0;
//     for (auto& pop : this->ions)
//     {
//         EXPECT_EQ(pop.domainParticles().size(), 100 * nbrPartPerCell);
//         EXPECT_EQ(pop.levelGhostParticles().size(), 2 * nbrPartPerCell);
//         EXPECT_EQ(pop.patchGhostParticles().size(), 0);
//         per_particle(pop.domainParticles(),
//                      [&](Particle<dim>& p) { deltaSum += PHARE::core::sum(p.delta()); });
//     }
//     EXPECT_DOUBLE_EQ(deltaSum, 103333.29999993632);

//     IonUpdater ionUpdater{};
//     this->update(ionUpdater, UpdaterMode::all);

//     std::size_t icellSum = 0, pCount = 0;
//     deltaSum          = 0;
//     double densitySum = 0, fluxXSum = 0, fluxYSum = 0, fluxZSum = 0;
//     for (auto& pop : this->ions)
//     {
//         densitySum += PHARE::core::sum(reduce(pop.density()));
//         fluxXSum += PHARE::core::sum(reduce(pop.flux()[0]));
//         fluxYSum += PHARE::core::sum(reduce(pop.flux()[1]));
//         fluxZSum += PHARE::core::sum(reduce(pop.flux()[2]));
//         per_particle(pop.domainParticles(),
//                      [&](Particle<dim>& p) { deltaSum += PHARE::core::sum(p.delta()); });
//         per_particle(pop.domainParticles(),
//                      [&](Particle<dim>& p) { icellSum += PHARE::core::sum(p.iCell()); });
//         pCount += pop.domainParticles().size();
//     }

//     EXPECT_EQ(pCount, 200000);
//     EXPECT_EQ(icellSum, 9900000);
//     EXPECT_DOUBLE_EQ(deltaSum, 123333.30000002828);
//     EXPECT_DOUBLE_EQ(densitySum, expectations[interp - 1]);
//     EXPECT_DOUBLE_EQ(fluxXSum, expectations[interp - 1]);
//     EXPECT_DOUBLE_EQ(fluxYSum, expectations[interp - 1]);
//     EXPECT_DOUBLE_EQ(fluxZSum, expectations[interp - 1]);
// }



} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
