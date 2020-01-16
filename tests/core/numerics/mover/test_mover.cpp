#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <iterator>

#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/ions/ion_population/particle_pack.h"
#include "core/data/ions/particle_initializers/particle_initializer_factory.h"

#include "core/numerics/ion_updater/ion_updater.h"

#include "simulator/phare_types.h"

//#include "initializer/data_provider.h"

using namespace PHARE::core;

/**

  - we need to create an ion object
        - with a dictionnay having every ion arguments to load particles
        - then use ion.setBuffer() to set the density ptr to a field ptr of the same type
        - then get a ref to the bulk velocity and use .setBuffer() on that to set vx,y,z
        - then for each ion population:
            - use setBuffer(name, field) to set the density
            - use setBuffer(name, pack) to set the particlePack
            - get the Flux VecField and use setBuffer() on that to set fx,y,z

  all fields (rhos) and vecfield components (v/f x,y,z) and particle packs
  need to be manually allocated.


  - create an electromag object
        - with a dictionnary with initialization functions
        - get E and B VecField ans use setBuffer() with manually allocated fields

  - create a layout consistent with allocated sizes


  - initialize electromag object
  - for each ion pop create a particle initializer and load the particles.

  - loading particles only loads domain particles
        - as is, some of the domain nodes will be incomplete and first ghost node too
        -  the deposit of levelGhostParticles and patchGhostParticle needs to be checked
        - so we assume levelGhost sides and patchGhost
  - check the moments are equal to prescribed values
  - save moments and particles
  - apply tests

 */



double density(double x)
{
    return 1.;
}

double vx(double /*x*/)
{
    return 0.;
}


double vy(double /*x*/)
{
    return 0.;
}


double vz(double /*x*/)
{
    return 0.;
}


double vthx(double /*x*/)
{
    return 0.1;
}


double vthy(double /*x*/)
{
    return 0.1;
}


double vthz(double /*x*/)
{
    return 0.1;
}


double bx(double x)
{
    (void)x;
    return 0.;
}

double by(double x)
{
    (void)x;
    return 0.;
}

double bz(double x)
{
    (void)x;
    return 1.;
}

double ex(double x)
{
    (void)x;
    return 0.;
}

double ey(double x)
{
    (void)x;
    return 1.;
}

double ez(double x)
{
    (void)x;
    return 0.;
}



using ScalarFunctionT = PHARE::initializer::ScalarFunction<1>;

PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;

    dict["simulation"]["solverPPC"]["pusher"]["name"] = std::string{"modified_boris"};

    dict["ions"]["name"]                                    = std::string{"ions"};
    dict["ions"]["nbrPopulations"]                          = int{2};
    dict["ions"]["pop0"]["name"]                            = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]                            = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["name"]    = std::string{"maxwellian"};
    dict["ions"]["pop0"]["particle_initializer"]["density"] = static_cast<ScalarFunctionT>(density);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<ScalarFunctionT>(vx);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<ScalarFunctionT>(vy);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<ScalarFunctionT>(vz);


    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<ScalarFunctionT>(vthx);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<ScalarFunctionT>(vthy);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<ScalarFunctionT>(vthz);


    dict["ions"]["pop0"]["particle_initializer"]["nbr_part_per_cell"] = int{100};
    dict["ions"]["pop0"]["particle_initializer"]["charge"]            = -1.;
    dict["ions"]["pop0"]["particle_initializer"]["basis"]             = std::string{"cartesian"};

    dict["ions"]["pop1"]["name"]                            = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"]                            = 1.;
    dict["ions"]["pop1"]["particle_initializer"]["name"]    = std::string{"maxwellian"};
    dict["ions"]["pop1"]["particle_initializer"]["density"] = static_cast<ScalarFunctionT>(density);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<ScalarFunctionT>(vx);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<ScalarFunctionT>(vy);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<ScalarFunctionT>(vz);


    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<ScalarFunctionT>(vthx);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<ScalarFunctionT>(vthy);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<ScalarFunctionT>(vthz);


    dict["ions"]["pop1"]["particle_initializer"]["nbr_part_per_cell"] = int{100};
    dict["ions"]["pop1"]["particle_initializer"]["charge"]            = -1.;
    dict["ions"]["pop1"]["particle_initializer"]["basis"]             = std::string{"cartesian"};

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["electric"]["initializer"]["x_component"] = static_cast<ScalarFunctionT>(ex);
    dict["electromag"]["electric"]["initializer"]["y_component"] = static_cast<ScalarFunctionT>(ey);
    dict["electromag"]["electric"]["initializer"]["z_component"] = static_cast<ScalarFunctionT>(ez);

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<ScalarFunctionT>(bx);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<ScalarFunctionT>(by);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<ScalarFunctionT>(bz);

    return dict;
}



template<std::size_t dim, std::size_t interporder>
struct DimInterp
{
    static constexpr auto dimension    = dim;
    static constexpr auto interp_order = interporder;
};



// the Electromag and Ions used in this test
// need their resources pointers (Fields and ParticleArrays) to set manually
// to buffers. ElectromagBuffer and IonsBuffer encapsulate these buffers



template<std::size_t dim, std::size_t interp_order>
struct ElectromagBuffers
{
    using PHARETypes = PHARE::PHARE_Types<dim, interp_order>;
    using Field      = typename PHARETypes::Field_t;
    using GridLayout = typename PHARETypes::GridLayout_t;
    using Electromag = typename PHARETypes::Electromag_t;

    Field Bx, By, Bz;
    Field Ex, Ey, Ez;

    ElectromagBuffers(GridLayout const& layout)
        : Bx{"EM_B_x", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"EM_B_y", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"EM_B_z", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"EM_E_x", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Ey{"EM_E_x", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Ez{"EM_E_x", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
    {
    }

    ElectromagBuffers(ElectromagBuffers const& source, GridLayout const& layout)
        : Bx{"EM_B_x", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"EM_B_y", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"EM_B_z", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"EM_E_x", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Ey{"EM_E_x", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Ez{"EM_E_x", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
    {
        Bx.copyData(source.Bx);
        By.copyData(source.By);
        Bz.copyData(source.Bz);

        Ex.copyData(source.Ex);
        Ey.copyData(source.Ey);
        Ez.copyData(source.Ez);
    }


    void setBuffers(Electromag& EM)
    {
        EM.B.setBuffer("EM_B_x", &Bx);
        EM.B.setBuffer("EM_B_y", &By);
        EM.B.setBuffer("EM_B_z", &Bz);

        EM.E.setBuffer("EM_E_x", &Ex);
        EM.E.setBuffer("EM_E_y", &Ey);
        EM.E.setBuffer("EM_E_z", &Ez);
    }
};




template<std::size_t dim, std::size_t interp_order>
struct IonsBuffers
{
    using PHARETypes                 = PHARE::PHARE_Types<dim, interp_order>;
    using Field                      = typename PHARETypes::Field_t;
    using GridLayout                 = typename PHARETypes::GridLayout_t;
    using Ions                       = typename PHARETypes::Ions_t;
    using ParticleArray              = typename PHARETypes::ParticleArray_t;
    using ParticleInitializerFactory = typename PHARETypes::ParticleInitializerFactory;

    Field ionDensity;
    Field protonDensity;
    Field alphaDensity;
    Field protonFx;
    Field protonFy;
    Field protonFz;
    Field alphaFx;
    Field alphaFy;
    Field alphaFz;
    Field Vx;
    Field Vy;
    Field Vz;

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
        : ionDensity{"ions_rho", HybridQuantity::Scalar::rho,
                     layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonDensity{"ions_protons_rho", HybridQuantity::Scalar::rho,
                        layout.allocSize(HybridQuantity::Scalar::rho)}
        , alphaDensity{"ions_alpha_rho", HybridQuantity::Scalar::rho,
                       layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonFx{"ions_protons_flux_x", HybridQuantity::Scalar::Vx,
                   layout.allocSize(HybridQuantity::Scalar::Vx)}
        , protonFy{"ions_protons_flux_y", HybridQuantity::Scalar::Vy,
                   layout.allocSize(HybridQuantity::Scalar::Vy)}
        , protonFz{"ions_protons_flux_z", HybridQuantity::Scalar::Vz,
                   layout.allocSize(HybridQuantity::Scalar::Vz)}
        , alphaFx{"ions_alpha_flux_x", HybridQuantity::Scalar::Vx,
                  layout.allocSize(HybridQuantity::Scalar::Vx)}
        , alphaFy{"ions_alpha_flux_y", HybridQuantity::Scalar::Vy,
                  layout.allocSize(HybridQuantity::Scalar::Vy)}
        , alphaFz{"ions_alpha_flux_z", HybridQuantity::Scalar::Vz,
                  layout.allocSize(HybridQuantity::Scalar::Vz)}
        , Vx{"ions_bulkVel_x", HybridQuantity::Scalar::Vx,
             layout.allocSize(HybridQuantity::Scalar::Vx)}
        , Vy{"ions_bulkVel_y", HybridQuantity::Scalar::Vy,
             layout.allocSize(HybridQuantity::Scalar::Vy)}
        , Vz{"ions_bulkVel_z", HybridQuantity::Scalar::Vz,
             layout.allocSize(HybridQuantity::Scalar::Vz)}

    {
        protonPack.domainParticles        = &protonDomain;
        protonPack.patchGhostParticles    = &protonPatchGhost;
        protonPack.levelGhostParticles    = &protonLevelGhost;
        protonPack.levelGhostParticlesOld = &protonLevelGhostOld;
        protonPack.levelGhostParticlesNew = &protonLevelGhostNew;
        alphaPack.domainParticles         = &alphaDomain;
        alphaPack.patchGhostParticles     = &alphaPatchGhost;
        alphaPack.levelGhostParticles     = &alphaLevelGhost;
        alphaPack.levelGhostParticlesOld  = &alphaLevelGhostOld;
        alphaPack.levelGhostParticlesNew  = &alphaLevelGhostNew;
    }


    IonsBuffers(IonsBuffers const& source, GridLayout const& layout)
        : ionDensity{"ions_rho", HybridQuantity::Scalar::rho,
                     layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonDensity{"ions_protons_rho", HybridQuantity::Scalar::rho,
                        layout.allocSize(HybridQuantity::Scalar::rho)}
        , alphaDensity{"ions_alpha_rho", HybridQuantity::Scalar::rho,
                       layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonFx{"ions_protons_flux_x", HybridQuantity::Scalar::Vx,
                   layout.allocSize(HybridQuantity::Scalar::Vx)}
        , protonFy{"ions_protons_flux_y", HybridQuantity::Scalar::Vy,
                   layout.allocSize(HybridQuantity::Scalar::Vy)}
        , protonFz{"ions_protons_flux_z", HybridQuantity::Scalar::Vz,
                   layout.allocSize(HybridQuantity::Scalar::Vz)}
        , alphaFx{"ions_alpha_flux_x", HybridQuantity::Scalar::Vx,
                  layout.allocSize(HybridQuantity::Scalar::Vx)}
        , alphaFy{"ions_alpha_flux_y", HybridQuantity::Scalar::Vy,
                  layout.allocSize(HybridQuantity::Scalar::Vy)}
        , alphaFz{"ions_alpha_flux_z", HybridQuantity::Scalar::Vz,
                  layout.allocSize(HybridQuantity::Scalar::Vz)}
        , Vx{"ions_bulkVel_x", HybridQuantity::Scalar::Vx,
             layout.allocSize(HybridQuantity::Scalar::Vx)}
        , Vy{"ions_bulkVel_y", HybridQuantity::Scalar::Vy,
             layout.allocSize(HybridQuantity::Scalar::Vy)}
        , Vz{"ions_bulkVel_z", HybridQuantity::Scalar::Vz,
             layout.allocSize(HybridQuantity::Scalar::Vz)}
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

    {
        ionDensity.copyData(source.ionDensity);
        protonDensity.copyData(source.protonDensity);
        alphaDensity.copyData(source.alphaDensity);

        protonFx.copyData(source.protonFx);
        protonFy.copyData(source.protonFy);
        protonFz.copyData(source.protonFz);

        alphaFx.copyData(source.alphaFx);
        alphaFy.copyData(source.alphaFy);
        alphaFz.copyData(source.alphaFz);

        Vx.copyData(source.Vx);
        Vy.copyData(source.Vy);
        Vz.copyData(source.Vz);

        protonPack.domainParticles        = &protonDomain;
        protonPack.patchGhostParticles    = &protonPatchGhost;
        protonPack.levelGhostParticles    = &protonLevelGhost;
        protonPack.levelGhostParticlesOld = &protonLevelGhostOld;
        protonPack.levelGhostParticlesNew = &protonLevelGhostNew;
        alphaPack.domainParticles         = &alphaDomain;
        alphaPack.patchGhostParticles     = &alphaPatchGhost;
        alphaPack.levelGhostParticles     = &alphaLevelGhost;
        alphaPack.levelGhostParticlesOld  = &alphaLevelGhostOld;
        alphaPack.levelGhostParticlesNew  = &alphaLevelGhostNew;
    }

    void setBuffers(Ions& ions)
    {
        ions.setBuffer("ions_rho", &ionDensity);
        auto& v = ions.velocity();
        v.setBuffer("ions_bulkVel_x", &Vx);
        v.setBuffer("ions_bulkVel_y", &Vy);
        v.setBuffer("ions_bulkVel_z", &Vz);

        auto& populations = ions.getRunTimeResourcesUserList();

        populations[0].setBuffer("ions_protons_rho", &protonDensity);
        populations[0].flux().setBuffer("ions_protons_flux_x", &protonFx);
        populations[0].flux().setBuffer("ions_protons_flux_y", &protonFy);
        populations[0].flux().setBuffer("ions_protons_flux_z", &protonFz);


        populations[0].setBuffer("ions_protons", &protonPack);

        populations[1].setBuffer("ions_alpha_rho", &alphaDensity);
        populations[1].flux().setBuffer("ions_alpha_flux_x", &alphaFx);
        populations[1].flux().setBuffer("ions_alpha_flux_y", &alphaFy);
        populations[1].flux().setBuffer("ions_alpha_flux_z", &alphaFz);


        populations[1].setBuffer("ions_alpha", &alphaPack);
    }
};




template<typename DimInterpT>
struct IonUpdaterTest : public ::testing::Test
{
    static constexpr auto dim          = DimInterpT::dimension;
    static constexpr auto interp_order = DimInterpT::interp_order;
    using PHARETypes                   = PHARE::PHARE_Types<dim, interp_order>;
    using Ions                         = typename PHARETypes::Ions_t;
    using Electromag                   = typename PHARETypes::Electromag_t;
    using ParticleArray                = typename PHARETypes::ParticleArray_t;
    using ParticleInitializerFactory   = typename PHARETypes::ParticleInitializerFactory;


    // grid configuration
    std::array<int, dim> ncells;
    GridLayout<GridLayoutImplYee<dim, interp_order>> layout;


    // data for electromagnetic fields
    using Field    = typename PHARETypes::Field_t;
    using VecField = typename PHARETypes::VecField_t;

    ElectromagBuffers<dim, interp_order> emBuffers;
    IonsBuffers<dim, interp_order> ionsBuffers;

    Electromag EM{createDict()["electromag"]};
    Ions ions{createDict()["ions"]};



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
            auto info                = pop.particleInitializerInfo();
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


                std::copy_if(std::begin(domainPart), std::end(domainPart),
                             std::back_inserter(levelGhostPartOld),
                             [&firstAMRCell](auto const& particle) {
                                 if constexpr (interp_order == 1)
                                 {
                                     return particle.iCell[0] == firstAMRCell[0];
                                 }
                                 else if constexpr (interp_order == 2 or interp_order == 3)
                                 {
                                     return (particle.iCell[0] == firstAMRCell[0])
                                            || (particle.iCell[0] == firstAMRCell[0] + 1);
                                 }
                             });

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
                std::transform(std::begin(levelGhostPartOld), std::end(levelGhostPartOld),
                               std::begin(levelGhostPartOld), [](auto& part) {
                                   if constexpr (interp_order == 2 or interp_order == 3)
                                   {
                                       part.iCell[0] = part.iCell[0] - 2;
                                   }
                                   else if constexpr (interp_order == 1)
                                   {
                                       part.iCell[0] = part.iCell[0] - 1;
                                   }
                                   return part;
                               });


                std::copy(std::begin(levelGhostPartOld), std::end(levelGhostPartOld),
                          std::back_inserter(levelGhostPartNew));


                std::copy(std::begin(levelGhostPartOld), std::end(levelGhostPartOld),
                          std::back_inserter(levelGhostPart));


                // now let's create patchGhostParticles on the right of the domain
                // by copying those on the last cell

                std::copy_if(std::begin(domainPart), std::end(domainPart),
                             std::back_inserter(patchGhostPart),
                             [&lastAMRCell](auto const& particle) {
                                 if constexpr (interp_order == 1)
                                 {
                                     return particle.iCell[0] == lastAMRCell[0];
                                 }
                                 else if constexpr (interp_order == 2 or interp_order == 3)
                                 {
                                     return (particle.iCell[0] == lastAMRCell[0])
                                            || (particle.iCell[0] == lastAMRCell[0] - 1);
                                 }
                             });


                std::transform(std::begin(patchGhostPart), std::end(patchGhostPart),
                               std::begin(patchGhostPart), [](auto& part) {
                                   if constexpr (interp_order == 2 or interp_order == 3)
                                   {
                                       part.iCell[0] = part.iCell[0] + 2;
                                   }
                                   else if constexpr (interp_order == 1)
                                   {
                                       part.iCell[0] = part.iCell[0] + 1;
                                   }
                                   return part;
                               });


            } // end 1D
        }     // end pop loop
    }
};



using DimInterps = ::testing::Types<DimInterp<1, 1>, DimInterp<1, 2>, DimInterp<1, 3>>;


TYPED_TEST_SUITE(IonUpdaterTest, DimInterps);


TYPED_TEST(IonUpdaterTest, momentsAreUpdatedButParticlesUnTouchedInMomentOnlyMode)
{
    if constexpr (TypeParam::dimension == 1)
    {
        IonUpdater ionUpdater{};

        ElectromagBuffers emBufferCpy{this->emBuffers, this->layout};
        IonsBuffers ionsBufferCpy{this->ionsBuffers, this->layout};

        ionUpdater.update(this->ions, this->EM, this->layout, UpdaterMode::moments_only);

        for (auto& pop : this->ions)
        {
            ASSERT_EQ(this->layout.nbrCells()[0] * 100, pop.domainParticles().size());
            ASSERT_EQ(100, pop.patchGhostParticles().size());
            ASSERT_EQ(100, pop.levelGhostParticlesOld().size());
            ASSERT_EQ(100, pop.levelGhostParticlesNew().size());
            ASSERT_EQ(100, pop.levelGhostParticles().size());
        }
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
