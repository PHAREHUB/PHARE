#include "gtest/gtest.h"
#include <core/data/grid/gridlayout.hpp>

#include "phare_core.hpp"

#include "core/numerics/ion_updater/ion_updater.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <numeric>
#include <algorithm>


using namespace PHARE::core;


using Param  = PHARE::core::Span<double const> const&;
using Return = std::shared_ptr<Span<double>>;

Return density(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), 1);
}

Return vx(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), .1);
}

Return vy(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), .1);
}

Return vz(Param x)
{
    return std::make_shared<VectorSpan<double>>(x.size(), .1);
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




int nbrPartPerCell = 1000;

using InitFunctionT = PHARE::initializer::InitFunction<1>;

PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;

    dict["simulation"]["algo"]["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};

    dict["ions"]["nbrPopulations"] = std::size_t{2};

    dict["ions"]["pop0"]["name"]           = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]           = 1.;
    auto& pop0_initializer                 = dict["ions"]["pop0"]["particle_initializer"];
    pop0_initializer["name"]               = std::string{"maxwellian"};
    pop0_initializer["density"]            = static_cast<InitFunctionT>(density);
    pop0_initializer["bulk_velocity_x"]    = static_cast<InitFunctionT>(vx);
    pop0_initializer["bulk_velocity_y"]    = static_cast<InitFunctionT>(vy);
    pop0_initializer["bulk_velocity_z"]    = static_cast<InitFunctionT>(vz);
    pop0_initializer["thermal_velocity_x"] = static_cast<InitFunctionT>(vthx);
    pop0_initializer["thermal_velocity_y"] = static_cast<InitFunctionT>(vthy);
    pop0_initializer["thermal_velocity_z"] = static_cast<InitFunctionT>(vthz);
    pop0_initializer["nbr_part_per_cell"]  = int{nbrPartPerCell};
    pop0_initializer["charge"]             = 1.;
    pop0_initializer["basis"]              = std::string{"cartesian"};

    dict["ions"]["pop1"]["name"]           = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"]           = 1.;
    auto& pop1_initializer                 = dict["ions"]["pop1"]["particle_initializer"];
    pop1_initializer["name"]               = std::string{"maxwellian"};
    pop1_initializer["density"]            = static_cast<InitFunctionT>(density);
    pop1_initializer["bulk_velocity_x"]    = static_cast<InitFunctionT>(vx);
    pop1_initializer["bulk_velocity_y"]    = static_cast<InitFunctionT>(vy);
    pop1_initializer["bulk_velocity_z"]    = static_cast<InitFunctionT>(vz);
    pop1_initializer["thermal_velocity_x"] = static_cast<InitFunctionT>(vthx);
    pop1_initializer["thermal_velocity_y"] = static_cast<InitFunctionT>(vthy);
    pop1_initializer["thermal_velocity_z"] = static_cast<InitFunctionT>(vthz);
    pop1_initializer["nbr_part_per_cell"]  = int{nbrPartPerCell};
    pop1_initializer["charge"]             = 1.;
    pop1_initializer["basis"]              = std::string{"cartesian"};

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<InitFunctionT>(bx);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<InitFunctionT>(by);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<InitFunctionT>(bz);

    return dict;
}
static auto init_dict = createDict();




// the Electromag and Ions used in this test
// need their resources pointers (Fields and ParticleArrays) to set manually
// to buffers. ElectromagBuffer and IonsBuffer encapsulate these buffers



template<std::size_t dim, std::size_t interp_order>
struct ElectromagBuffers
{
    constexpr static PHARE::SimOpts opts{dim, interp_order};
    using PHARETypes       = PHARE::core::PHARE_Types<opts>;
    using Grid             = PHARETypes::Hybrid::Grid_t;
    using GridLayout       = PHARETypes::Hybrid::GridLayout_t;
    using Electromag       = PHARETypes::Hybrid::Electromag_t;
    using UsableVecFieldND = UsableVecField<dim>;

    UsableVecFieldND B, E;

    ElectromagBuffers(GridLayout const& layout)
        : B{"EM_B", layout, HybridQuantity::Vector::B}
        , E{"EM_E", layout, HybridQuantity::Vector::E}
    {
    }

    ElectromagBuffers(ElectromagBuffers const& source, GridLayout const& layout)
        : ElectromagBuffers{layout}
    {
        B.copyData(source.B);
        E.copyData(source.E);
    }


    void setBuffers(Electromag& EM)
    {
        B.set_on(EM.B);
        E.set_on(EM.E);
    }
};




template<std::size_t dim, std::size_t interp_order>
struct IonsBuffers
{
    constexpr static PHARE::SimOpts opts{dim, interp_order};
    using PHARETypes                 = PHARE_Types<opts>;
    using UsableVecFieldND           = UsableVecField<dim>;
    using Grid                       = PHARETypes::Hybrid::Grid_t;
    using GridLayout                 = PHARETypes::Hybrid::GridLayout_t;
    using Ions                       = PHARETypes::Hybrid::Ions_t;
    using ParticleArray              = PHARETypes::Hybrid::ParticleArray_t;
    using ParticleInitializerFactory = PHARETypes::Hybrid::ParticleInitializerFactory_t;

    Grid ionChargeDensity;
    Grid ionMassDensity;
    Grid protonParticleDensity;
    Grid protonChargeDensity;
    Grid alphaParticleDensity;
    Grid alphaChargeDensity;

    UsableVecFieldND protonF, alphaF, Vi;
    UsableTensorField<dim> M, alpha_M, protons_M;

    static constexpr int ghostSafeMapLayer = ghostWidthForParticles<interp_order>() + 1;

    ParticleArray protonDomain;
    ParticleArray protonPatchGhost;
    ParticleArray protonLevelGhost;
    ParticleArray protonLevelGhostOld;
    ParticleArray protonLevelGhostNew;

    ParticleArray alphaDomain;
    ParticleArray alphaPatchGhost;
    ParticleArray alphaLevelGhost;
    ParticleArray alphaLevelGhostOld;
    ParticleArray alphaLevelGhostNew;

    ParticlesPack<ParticleArray> protonPack;
    ParticlesPack<ParticleArray> alphaPack;

    IonsBuffers(GridLayout const& layout)
        : ionChargeDensity{"chargeDensity", HybridQuantity::Scalar::rho,
                           layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , ionMassDensity{"massDensity", HybridQuantity::Scalar::rho,
                         layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , protonParticleDensity{"protons_particleDensity", HybridQuantity::Scalar::rho,
                                layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , protonChargeDensity{"protons_chargeDensity", HybridQuantity::Scalar::rho,
                              layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , alphaParticleDensity{"alpha_particleDensity", HybridQuantity::Scalar::rho,
                               layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , alphaChargeDensity{"alpha_chargeDensity", HybridQuantity::Scalar::rho,
                             layout.allocSize(HybridQuantity::Scalar::rho), 0.}
        , protonF{"protons_flux", layout, HybridQuantity::Vector::V}
        , alphaF{"alpha_flux", layout, HybridQuantity::Vector::V}
        , Vi{"bulkVel", layout, HybridQuantity::Vector::V}
        , M{"momentumTensor", layout, HybridQuantity::Tensor::M}
        , alpha_M{"alpha_momentumTensor", layout, HybridQuantity::Tensor::M}
        , protons_M{"protons_momentumTensor", layout, HybridQuantity::Tensor::M}
        , protonDomain{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonPatchGhost{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonLevelGhost{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonLevelGhostOld{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonLevelGhostNew{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaDomain{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaPatchGhost{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaLevelGhost{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaLevelGhostOld{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , alphaLevelGhostNew{grow(layout.AMRBox(), ghostSafeMapLayer)}
        , protonPack{"protons",         &protonDomain,        &protonPatchGhost,
                     &protonLevelGhost, &protonLevelGhostOld, &protonLevelGhostNew}
        , alphaPack{"alpha",          &alphaDomain,        &alphaPatchGhost,
                    &alphaLevelGhost, &alphaLevelGhostOld, &alphaLevelGhostNew}
    {
    }


    IonsBuffers(IonsBuffers const& source, GridLayout const& layout)
        : ionChargeDensity{"chargeDensity", HybridQuantity::Scalar::rho,
                           layout.allocSize(HybridQuantity::Scalar::rho)}
        , ionMassDensity{"massDensity", HybridQuantity::Scalar::rho,
                         layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonParticleDensity{"protons_particleDensity", HybridQuantity::Scalar::rho,
                                layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonChargeDensity{"protons_chargeDensity", HybridQuantity::Scalar::rho,
                              layout.allocSize(HybridQuantity::Scalar::rho)}
        , alphaParticleDensity{"alpha_particleDensity", HybridQuantity::Scalar::rho,
                               layout.allocSize(HybridQuantity::Scalar::rho)}
        , alphaChargeDensity{"alpha_chargeDensity", HybridQuantity::Scalar::rho,
                             layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonF{"protons_flux", layout, HybridQuantity::Vector::V}
        , alphaF{"alpha_flux", layout, HybridQuantity::Vector::V}
        , Vi{"bulkVel", layout, HybridQuantity::Vector::V}
        , M{"momentumTensor", layout, HybridQuantity::Tensor::M}
        , alpha_M{"alpha_momentumTensor", layout, HybridQuantity::Tensor::M}
        , protons_M{"protons_momentumTensor", layout, HybridQuantity::Tensor::M}
        , protonDomain{source.protonDomain}
        , protonPatchGhost{source.protonPatchGhost}
        , protonLevelGhost{source.protonLevelGhost}
        , protonLevelGhostOld{source.protonLevelGhostOld}
        , protonLevelGhostNew{source.protonLevelGhostNew}
        , alphaDomain{source.alphaDomain}
        , alphaPatchGhost{source.alphaPatchGhost}
        , alphaLevelGhost{source.alphaLevelGhost}
        , alphaLevelGhostOld{source.alphaLevelGhostOld}
        , alphaLevelGhostNew{source.alphaLevelGhostNew}
        , protonPack{"protons",         &protonDomain,        &protonPatchGhost,
                     &protonLevelGhost, &protonLevelGhostOld, &protonLevelGhostNew}
        , alphaPack{"alpha",          &alphaDomain,        &alphaPatchGhost,
                    &alphaLevelGhost, &alphaLevelGhostOld, &alphaLevelGhostNew}

    {
        ionChargeDensity.copyData(source.ionChargeDensity);
        ionMassDensity.copyData(source.ionMassDensity);
        protonParticleDensity.copyData(source.protonParticleDensity);
        protonChargeDensity.copyData(source.protonChargeDensity);
        alphaParticleDensity.copyData(source.alphaParticleDensity);
        alphaChargeDensity.copyData(source.alphaChargeDensity);

        protonF.copyData(source.protonF);
        alphaF.copyData(source.alphaF);
        Vi.copyData(source.Vi);
    }

    void setBuffers(Ions& ions)
    {
        {
            auto const& [V, m, cd, md] = ions.getCompileTimeResourcesViewList();
            Vi.set_on(V);
            M.set_on(m);
            cd.setBuffer(&ionChargeDensity);
            md.setBuffer(&ionMassDensity);
        }

        auto& pops = ions.getRunTimeResourcesViewList();
        {
            auto const& [F, M, d, c, particles] = pops[0].getCompileTimeResourcesViewList();
            d.setBuffer(&protonParticleDensity);
            c.setBuffer(&protonChargeDensity);
            protons_M.set_on(M);
            protonF.set_on(F);
            particles.setBuffer(&protonPack);
        }

        {
            auto const& [F, M, d, c, particles] = pops[1].getCompileTimeResourcesViewList();
            d.setBuffer(&alphaParticleDensity);
            c.setBuffer(&alphaChargeDensity);
            alpha_M.set_on(M);
            alphaF.set_on(F);
            particles.setBuffer(&alphaPack);
        }
    }
};


template<auto opts_, std::uint16_t impl_>
struct TestParam
{
    auto static constexpr opts         = opts_;
    auto static constexpr impl         = impl_;
    static constexpr auto dimension    = opts.dimension;
    static constexpr auto interp_order = opts.interp_order;
};


template<auto opts, std::size_t impl>
struct updater_test_bits;

template<auto opts>
struct updater_test_bits<opts, 0>
{
    using PHARETypes = PHARE::core::PHARE_Types<opts>;
    using Ions       = PHARETypes::Hybrid::Ions_t;
    using IonUpdater = IonUpdater0<Ions>;
};

template<auto opts>
struct updater_test_bits<opts, 1>
{
    using PHARETypes = PHARE::core::PHARE_Types<opts>;
    using Ions       = PHARETypes::Hybrid::Ions_t;
    using IonUpdater = IonUpdater1<Ions>;
};




template<typename TestParam_t>
struct IonUpdaterTest : public ::testing::Test,
                        public updater_test_bits<TestParam_t::opts, TestParam_t::impl>
{
    auto static constexpr opts         = TestParam_t::opts;
    static constexpr auto dim          = opts.dimension;
    static constexpr auto interp_order = opts.interp_order;
    using PHARETypes                   = PHARE::core::PHARE_Types<opts>;
    using Ions                         = PHARETypes::Hybrid::Ions_t;
    using Electromag                   = PHARETypes::Hybrid::Electromag_t;
    using ParticleArray                = PHARETypes::Hybrid::ParticleArray_t;
    using ParticleInitializerFactory   = PHARETypes::Hybrid::ParticleInitializerFactory_t;
    using GridLayout                   = PHARETypes::Hybrid::GridLayout_t;

    using basics       = updater_test_bits<TestParam_t::opts, TestParam_t::impl>;
    using IonUpdater_t = basics::IonUpdater;
    using Boxing_t     = UpdaterSelectionBoxing<GridLayout, ParticleArray>;


    double dt{0.01};

    // grid configuration
    std::array<int, dim> ncells;
    GridLayout layout;
    // assumes only level ghost cells
    Box<int, dim> const domainBox = layout.AMRBox();
    Box<int, dim> const ghostBox  = grow(domainBox, GridLayout::options.particle_ghost_width);
    Boxing_t const boxing{layout, remove(ghostBox, domainBox)};


    // data for electromagnetic fields
    using Field            = PHARETypes::Hybrid::Grid_t;
    using VecField         = PHARETypes::Hybrid::VecField_t;
    using UsableVecFieldND = UsableVecField<dim>;

    ElectromagBuffers<dim, interp_order> emBuffers;
    IonsBuffers<dim, interp_order> ionsBuffers;

    Electromag EM{init_dict["electromag"]};
    Ions ions{init_dict["ions"]};



    IonUpdaterTest()
        : ncells{100}
        , layout{{0.1}, {100u}, {{0.}}}
        , emBuffers{layout}
        , ionsBuffers{layout}
    {
        emBuffers.setBuffers(EM);
        ionsBuffers.setBuffers(ions);


        // ok all resources pointers are set to buffers
        // now let's initialize Electromag fields to user input functions
        // and ion population particles to user supplied moments


        EM.initialize(layout);
        for (auto& pop : ions)
        {
            auto const& info         = pop.particleInitializerInfo();
            auto particleInitializer = ParticleInitializerFactory::create(info);
            particleInitializer->loadParticles(pop.domainParticles(), layout);
        }


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
                for (auto const& part : domainPart)
                {
                    if constexpr (interp_order == 2 or interp_order == 3)
                    {
                        if (part.iCell[0] == firstAMRCell[0]
                            or part.iCell[0] == firstAMRCell[0] + 1)
                        {
                            auto p{part};
                            p.iCell[0] -= 2;
                            levelGhostPartOld.push_back(p);
                        }
                    }
                    else if constexpr (interp_order == 1)
                    {
                        if (part.iCell[0] == firstAMRCell[0])
                        {
                            auto p{part};
                            p.iCell[0] -= 1;
                            levelGhostPartOld.push_back(p);
                        }
                    }
                }


                std::copy(std::begin(levelGhostPartOld), std::end(levelGhostPartOld),
                          std::back_inserter(levelGhostPartNew));


                std::copy(std::begin(levelGhostPartOld), std::end(levelGhostPartOld),
                          std::back_inserter(levelGhostPart));


                EXPECT_GT(pop.domainParticles().size(), 0ull);
                EXPECT_GT(levelGhostPartOld.size(), 0ull);
                EXPECT_EQ(patchGhostPart.size(), 0);

            } // end 1D
        } // end pop loop
        depositParticles(ions, layout, Interpolator<dim, interp_order>{}, DomainDeposit{});


        depositParticles(ions, layout, Interpolator<dim, interp_order>{}, LevelGhostDeposit{});


        ions.computeChargeDensity();
        ions.computeBulkVelocity();
    } // end Ctor



    void fillIonsMomentsGhosts()
    {
        using Interpolator = IonUpdater_t::Interpolator_t;
        Interpolator interpolate;

        for (auto& pop : this->ions)
        {
            double alpha = 0.5;
            interpolate(makeIndexRange(pop.levelGhostParticlesNew()), pop.particleDensity(),
                        pop.chargeDensity(), pop.flux(), layout,
                        /*coef = */ alpha);


            interpolate(makeIndexRange(pop.levelGhostParticlesOld()), pop.particleDensity(),
                        pop.chargeDensity(), pop.flux(), layout,
                        /*coef = */ (1. - alpha));
        }
    }



    void checkMomentsHaveEvolved(IonsBuffers<dim, interp_order> const& ionsBufferCpy)
    {
        auto& populations = this->ions.getRunTimeResourcesViewList();

        auto& protonParticleDensity = populations[0].particleDensity();
        auto& protonChargeDensity   = populations[0].chargeDensity();
        auto& protonFx              = populations[0].flux().getComponent(Component::X);
        auto& protonFy              = populations[0].flux().getComponent(Component::Y);
        auto& protonFz              = populations[0].flux().getComponent(Component::Z);

        auto& alphaParticleDensity = populations[1].particleDensity();
        auto& alphaChargeDensity   = populations[1].chargeDensity();
        auto& alphaFx              = populations[1].flux().getComponent(Component::X);
        auto& alphaFy              = populations[1].flux().getComponent(Component::Y);
        auto& alphaFz              = populations[1].flux().getComponent(Component::Z);

        auto ix0 = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
        auto ix1 = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

        auto nonZero = [&](auto const& field) {
            auto sum = 0.;
            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                sum += std::abs(field(ix));
            }
            EXPECT_GT(sum, 0.);
        };

        auto check = [&](auto const& newField, auto const& originalField) {
            nonZero(newField);
            nonZero(originalField);
            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                auto evolution = std::abs(newField(ix) - originalField(ix));
                //  should check that moments are still compatible with user inputs also
                EXPECT_TRUE(evolution > 0.0);
                if (evolution <= 0.0)
                    std::cout << "after update : " << newField(ix)
                              << " before update : " << originalField(ix)
                              << " evolution : " << evolution << " ix : " << ix << "\n";
            }
        };

        check(protonParticleDensity, ionsBufferCpy.protonParticleDensity);
        check(protonChargeDensity, ionsBufferCpy.protonChargeDensity);
        check(protonFx, ionsBufferCpy.protonF(Component::X));
        check(protonFy, ionsBufferCpy.protonF(Component::Y));
        check(protonFz, ionsBufferCpy.protonF(Component::Z));

        check(alphaParticleDensity, ionsBufferCpy.alphaParticleDensity);
        check(alphaChargeDensity, ionsBufferCpy.alphaChargeDensity);
        check(alphaFx, ionsBufferCpy.alphaF(Component::X));
        check(alphaFy, ionsBufferCpy.alphaF(Component::Y));
        check(alphaFz, ionsBufferCpy.alphaF(Component::Z));

        check(ions.velocity().getComponent(Component::X), ionsBufferCpy.Vi(Component::X));
        check(ions.velocity().getComponent(Component::Y), ionsBufferCpy.Vi(Component::Y));
        check(ions.velocity().getComponent(Component::Z), ionsBufferCpy.Vi(Component::Z));
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

            auto const span   = make_const_span(x);
            auto functionXPtr = function(span); // keep alive
            EXPECT_EQ(functionXPtr->size(), (ix1 - ix0));

            auto& functionX = *functionXPtr;

            for (std::size_t i = 0; i < functionX.size(); i++)
            {
                auto ix   = ixes[i];
                auto diff = std::abs(density(ix) - functionX[i]);

                EXPECT_GE(0.07, diff);

                if (diff >= 0.07)
                    std::cout << "actual : " << density(ix) << " prescribed : " << functionX[i]
                              << " diff : " << diff << " ix : " << ix << "\n";
            }
        };

        auto& populations           = this->ions.getRunTimeResourcesViewList();
        auto& protonParticleDensity = populations[0].particleDensity();
        auto& alphaParticleDensity  = populations[1].particleDensity();

        check(protonParticleDensity, density);
        check(alphaParticleDensity, density);
    }
};



// clang-format off
using Permutations = ::testing::Types<
    TestParam<PHARE::SimOpts{1, 1}, 0>
  , TestParam<PHARE::SimOpts{1, 2}, 0>
  , TestParam<PHARE::SimOpts{1, 3}, 0>
  , TestParam<PHARE::SimOpts{1, 1}, 1>
  , TestParam<PHARE::SimOpts{1, 2}, 1>
  , TestParam<PHARE::SimOpts{1, 3}, 1>
>;
// clang-format on
TYPED_TEST_SUITE(IonUpdaterTest, Permutations, );



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
            if constexpr (TypeParam::interp_order == 1)
            {
                for (auto const& part : pop.levelGhostParticles())
                {
                    EXPECT_EQ(firstAMRCell[0] - 1, part.iCell[0]);
                }
            }
            else if constexpr (TypeParam::interp_order == 2 or TypeParam::interp_order == 3)
            {
                typename IonUpdaterTest<TypeParam>::ParticleArray copy{pop.levelGhostParticles()};
                auto firstInOuterMostCell = std::partition(
                    std::begin(copy), std::end(copy), [&firstAMRCell](auto const& particle) {
                        return particle.iCell[0] == firstAMRCell[0] - 1;
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
    typename IonUpdaterTest<TypeParam>::IonUpdater_t ionUpdater{
        init_dict["simulation"]["algo"]["ion_updater"]};

    IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};

    ionUpdater.updatePopulations(UpdaterMode::copy, this->ions, this->EM, this->boxing, this->dt);

    this->fillIonsMomentsGhosts();

    this->ions.update();


    auto& populations = this->ions.getRunTimeResourcesViewList();

    auto checkIsUnTouched = [](auto const& original, auto const& cpy) {
        // no particles should have moved, so none should have left the domain
        EXPECT_EQ(cpy.size(), original.size());
        for (std::size_t iPart = 0; iPart < original.size(); ++iPart)
        {
            EXPECT_EQ(cpy[iPart].iCell[0], original[iPart].iCell[0]);
            EXPECT_DOUBLE_EQ(cpy[iPart].delta[0], original[iPart].delta[0]);

            for (std::size_t iDir = 0; iDir < 3; ++iDir)
            {
                EXPECT_DOUBLE_EQ(cpy[iPart].v[iDir], original[iPart].v[iDir]);
            }
        }
    };

    checkIsUnTouched(populations[0].patchGhostParticles(), ionsBufferCpy.protonPatchGhost);
    checkIsUnTouched(populations[0].levelGhostParticles(), ionsBufferCpy.protonLevelGhost);
    checkIsUnTouched(populations[0].levelGhostParticlesOld(), ionsBufferCpy.protonLevelGhostOld);
    checkIsUnTouched(populations[0].levelGhostParticlesNew(), ionsBufferCpy.protonLevelGhostNew);

    checkIsUnTouched(populations[1].patchGhostParticles(), ionsBufferCpy.alphaPatchGhost);
    checkIsUnTouched(populations[1].levelGhostParticles(), ionsBufferCpy.alphaLevelGhost);
    checkIsUnTouched(populations[1].levelGhostParticlesOld(), ionsBufferCpy.alphaLevelGhost);
    checkIsUnTouched(populations[1].levelGhostParticlesNew(), ionsBufferCpy.alphaLevelGhost);
}




// TYPED_TEST(IonUpdaterTest, particlesAreChangedInParticlesAndMomentsMode)
//{
//    typename IonUpdaterTest<TypeParam>::IonUpdater_t
//    ionUpdater{init_dict["simulation"]["pusher"]};
//
//    IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};
//
//    ionUpdater.updatePopulations(this->ions, this->EM, this->boxing, this->dt,
//                                 UpdaterMode::particles_and_moments);
//
//    this->fillIonsMomentsGhosts();
//
//    this->ions.update();
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
    typename IonUpdaterTest<TypeParam>::IonUpdater_t ionUpdater{
        init_dict["simulation"]["algo"]["ion_updater"]};

    IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};

    ionUpdater.updatePopulations(UpdaterMode::ref, this->ions, this->EM, this->boxing, this->dt);

    this->fillIonsMomentsGhosts();

    this->ions.update();

    this->checkMomentsHaveEvolved(ionsBufferCpy);
    this->checkDensityIsAsPrescribed();
}




TYPED_TEST(IonUpdaterTest, momentsAreChangedInMomentsOnlyMode)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater_t ionUpdater{
        init_dict["simulation"]["algo"]["ion_updater"]};

    IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};

    ionUpdater.updatePopulations(UpdaterMode::copy, this->ions, this->EM, this->boxing, this->dt);

    this->fillIonsMomentsGhosts();

    this->ions.update();

    this->checkMomentsHaveEvolved(ionsBufferCpy);
    this->checkDensityIsAsPrescribed();
}



TYPED_TEST(IonUpdaterTest, thatNoNaNsExistOnPhysicalNodesMoments)
{
    typename IonUpdaterTest<TypeParam>::IonUpdater_t ionUpdater{
        init_dict["simulation"]["algo"]["ion_updater"]};

    ionUpdater.updatePopulations(UpdaterMode::copy, this->ions, this->EM, this->boxing, this->dt);

    this->fillIonsMomentsGhosts();

    this->ions.update();

    auto ix0 = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto ix1 = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    for (auto& pop : this->ions)
    {
        for (auto ix = ix0; ix <= ix1; ++ix)
        {
            auto& density = pop.particleDensity();
            auto& flux    = pop.flux();

            auto& fx = flux.getComponent(Component::X);
            auto& fy = flux.getComponent(Component::Y);
            auto& fz = flux.getComponent(Component::Z);

            EXPECT_FALSE(std::isnan(density(ix)));
            EXPECT_FALSE(std::isnan(fx(ix)));
            EXPECT_FALSE(std::isnan(fy(ix)));
            EXPECT_FALSE(std::isnan(fz(ix)));
        }
    }
}




// Overlapping same-level patches: cells in `nonOwnedBoxes` belong to another patch.
// Particles move only along x, at a fixed fraction of a cell per step (E = B = 0).
template<typename TestParam_t>
struct IonUpdaterOwnershipTest : public IonUpdaterTest<TestParam_t>
{
    using Super    = IonUpdaterTest<TestParam_t>;
    using Box_t    = Box<int, Super::dim>;
    using Boxing_t = Super::Boxing_t;

    static constexpr double cellsPerStep = 0.5;

    IonUpdaterOwnershipTest()
    {
        this->EM.E.zero();
        this->EM.B.zero();
    }

    Boxing_t makeBoxing(std::vector<Box_t> const& nonOwned) const
    {
        return {this->layout,
                {grow(this->layout.AMRBox(), Super::GridLayout::options.particle_ghost_width)},
                nonOwned};
    }

    static Box_t box(int lower, int upper) { return Box_t{{lower}, {upper}}; }

    void setVelocityX(auto& particles, double const v)
    {
        auto const dx = this->layout.meshSize()[0];
        for (auto& particle : particles)
            particle.v[0] = v * cellsPerStep * dx / this->dt;
    }

    // what level initialization does on overlapping patches
    void eraseNonOwned(std::vector<Box_t> const& nonOwned)
    {
        for (auto& pop : this->ions)
        {
            auto& particles = pop.domainParticles();
            auto range      = makeIndexRange(particles);
            auto const kept = particles.partition(
                range, [&](auto const& cell) { return !isIn(Point{cell}, nonOwned); });
            particles.erase(makeRange(particles, kept.iend(), particles.size()));
        }
    }

    static std::size_t count(auto const& particles, auto&& predicate)
    {
        return sum_from(particles, [&](auto const& p) { return predicate(p) ? 1ul : 0ul; });
    }

    // particles of `cell` whose displacement takes them into the next cell
    static auto crossing(auto const& particles, int const cell)
    {
        return count(particles, [&](auto const& p) {
            return p.iCell[0] == cell and p.delta[0] + cellsPerStep >= 1.;
        });
    }

    static double weights(auto const& particles)
    {
        return sum_from(particles, [](auto const& particle) { return particle.weight; });
    }

    static double total(auto const& field) { return sum(field); }

    // level ghosts move into the first domain cell, which is owned or not
    void checkLevelGhostExport(bool const withNonOwned)
    {
        typename Super::IonUpdater_t ionUpdater{init_dict["simulation"]["algo"]["ion_updater"]};

        auto const nonOwned = withNonOwned ? std::vector{box(0, 9)} : std::vector<Box_t>{};
        auto const boxing   = this->makeBoxing(nonOwned);
        eraseNonOwned(nonOwned);

        std::vector<std::size_t> nbrDomain, nbrEntering;
        for (auto& pop : this->ions)
        {
            setVelocityX(pop.domainParticles(), 0);
            setVelocityX(pop.levelGhostParticles(), 1);
            nbrDomain.push_back(pop.domainParticles().size());
            nbrEntering.push_back(crossing(pop.levelGhostParticles(), -1));
            ASSERT_GT(nbrEntering.back(), 0u);
        }

        ionUpdater.updatePopulations(UpdaterMode::ref, this->ions, this->EM, boxing, this->dt);

        std::size_t ipop = 0;
        for (auto& pop : this->ions)
        {
            auto const& domain = pop.domainParticles();
            EXPECT_EQ(domain.size(), nbrDomain[ipop] + (withNonOwned ? 0 : nbrEntering[ipop]));
            EXPECT_EQ(count(domain, [&](auto const& p) { return !boxing.isOwned(p.iCell); }), 0u);
            ++ipop;
        }
    }

    void checkLevelGhostDeposit(bool const withNonOwned)
    {
        typename Super::IonUpdater_t ionUpdater{init_dict["simulation"]["algo"]["ion_updater"]};

        auto const nonOwned = withNonOwned ? std::vector{box(0, 9)} : std::vector<Box_t>{};
        auto const boxing   = this->makeBoxing(nonOwned);
        eraseNonOwned(nonOwned);

        std::vector<double> expected;
        for (auto& pop : this->ions)
        {
            setVelocityX(pop.domainParticles(), 0);
            setVelocityX(pop.levelGhostParticles(), 1);

            double enteringWeights = 0;
            for (auto const& p : pop.levelGhostParticles())
                if (p.iCell[0] == -1 and p.delta[0] + cellsPerStep >= 1.)
                    enteringWeights += p.weight;
            ASSERT_GT(enteringWeights, 0.);

            expected.push_back(weights(pop.domainParticles())
                               + (withNonOwned ? 0. : enteringWeights));
        }

        ionUpdater.updatePopulations(UpdaterMode::copy, this->ions, this->EM, boxing, this->dt);

        std::size_t ipop = 0;
        for (auto& pop : this->ions)
        {
            EXPECT_NEAR(total(pop.particleDensity()), expected[ipop], 1e-10 * expected[ipop]);
            ++ipop;
        }
    }
};

TYPED_TEST_SUITE(IonUpdaterOwnershipTest, Permutations, );



TYPED_TEST(IonUpdaterOwnershipTest, particleEnteringNonOwnedCellLeavesDomainAndIsDepositedOnce)
{
    typename TestFixture::IonUpdater_t ionUpdater{init_dict["simulation"]["algo"]["ion_updater"]};

    auto const nonOwned = std::vector{TestFixture::box(60, 99)};
    auto const boxing   = this->makeBoxing(nonOwned);
    this->eraseNonOwned(nonOwned);

    std::vector<std::size_t> nbrDomain, nbrCrossing;
    for (auto& pop : this->ions)
    {
        this->setVelocityX(pop.domainParticles(), 1);
        this->setVelocityX(pop.levelGhostParticles(), 0);
        nbrDomain.push_back(pop.domainParticles().size());
        nbrCrossing.push_back(TestFixture::crossing(pop.domainParticles(), 59));
        ASSERT_GT(nbrCrossing.back(), 0u);
    }

    ionUpdater.updatePopulations(UpdaterMode::ref, this->ions, this->EM, boxing, this->dt);

    std::size_t ipop = 0;
    for (auto& pop : this->ions)
    {
        auto const& domain     = pop.domainParticles();
        auto const& patchGhost = pop.patchGhostParticles();

        EXPECT_EQ(domain.size(), nbrDomain[ipop] - nbrCrossing[ipop]);
        EXPECT_EQ(
            TestFixture::count(domain, [&](auto const& p) { return !boxing.isOwned(p.iCell); }),
            0u);

        // the leaving particles wait in the patch ghost array for the exchange
        EXPECT_EQ(patchGhost.size(), nbrCrossing[ipop]);
        EXPECT_EQ(TestFixture::count(patchGhost, [](auto const& p) { return p.iCell[0] == 60; }),
                  nbrCrossing[ipop]);

        auto const deposited = TestFixture::total(pop.particleDensity());
        auto const expected  = TestFixture::weights(domain) + TestFixture::weights(patchGhost);
        EXPECT_NEAR(deposited, expected, 1e-10 * expected);
        ++ipop;
    }
}



TYPED_TEST(IonUpdaterOwnershipTest, levelGhostEnteringOwnedCellIsExported)
{
    this->checkLevelGhostExport(/*withNonOwned=*/false);
}

TYPED_TEST(IonUpdaterOwnershipTest, levelGhostEnteringNonOwnedCellIsNotExported)
{
    this->checkLevelGhostExport(/*withNonOwned=*/true);
}

TYPED_TEST(IonUpdaterOwnershipTest, levelGhostEnteringOwnedCellIsDepositedInMomentsOnlyMode)
{
    this->checkLevelGhostDeposit(/*withNonOwned=*/false);
}

TYPED_TEST(IonUpdaterOwnershipTest, levelGhostEnteringNonOwnedCellIsNotDepositedInMomentsOnlyMode)
{
    this->checkLevelGhostDeposit(/*withNonOwned=*/true);
}



TYPED_TEST(IonUpdaterOwnershipTest, ownershipFollowsDomainAndNonOwnedBoxes)
{
    auto const boxing = this->makeBoxing({TestFixture::box(10, 19), TestFixture::box(50, 50)});
    auto const owned  = [&](int i) { return boxing.isOwned(std::array{i}); };

    EXPECT_FALSE(owned(-1));
    EXPECT_TRUE(owned(0));
    EXPECT_TRUE(owned(9));
    EXPECT_FALSE(owned(10));
    EXPECT_FALSE(owned(19));
    EXPECT_TRUE(owned(20));
    EXPECT_FALSE(owned(50));
    EXPECT_TRUE(owned(99));
    EXPECT_FALSE(owned(100));
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
