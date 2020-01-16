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




template<typename DimInterpT>
struct IonMoverTest : public ::testing::Test
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

    Field Bx{"EM_B_x", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)};
    Field By{"EM_B_y", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)};
    Field Bz{"EM_B_z", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)};

    Field Ex{"EM_E_x", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)};
    Field Ey{"EM_E_x", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)};
    Field Ez{"EM_E_x", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)};

    Electromag EM{createDict()["electromag"]};


    // data for ions
    Ions ions{createDict()["ions"]};
    Field ionDensity{"ions_rho", HybridQuantity::Scalar::rho,
                     layout.allocSize(HybridQuantity::Scalar::rho)};

    Field protonDensity{"ions_protons_rho", HybridQuantity::Scalar::rho,
                        layout.allocSize(HybridQuantity::Scalar::rho)};

    Field alphaDensity{"ions_alpha_rho", HybridQuantity::Scalar::rho,
                       layout.allocSize(HybridQuantity::Scalar::rho)};

    Field protonFx{"ions_protons_flux_x", HybridQuantity::Scalar::Vx,
                   layout.allocSize(HybridQuantity::Scalar::Vx)};
    Field protonFy{"ions_protons_flux_y", HybridQuantity::Scalar::Vy,
                   layout.allocSize(HybridQuantity::Scalar::Vy)};
    Field protonFz{"ions_protons_flux_z", HybridQuantity::Scalar::Vz,
                   layout.allocSize(HybridQuantity::Scalar::Vz)};

    Field alphaFx{"ions_alpha_flux_x", HybridQuantity::Scalar::Vx,
                  layout.allocSize(HybridQuantity::Scalar::Vx)};
    Field alphaFy{"ions_alpha_flux_y", HybridQuantity::Scalar::Vy,
                  layout.allocSize(HybridQuantity::Scalar::Vy)};
    Field alphaFz{"ions_alpha_flux_z", HybridQuantity::Scalar::Vz,
                  layout.allocSize(HybridQuantity::Scalar::Vz)};

    Field Vx{"ions_bulkVel_x", HybridQuantity::Scalar::Vx,
             layout.allocSize(HybridQuantity::Scalar::Vx)};
    Field Vy{"ions_bulkVel_y", HybridQuantity::Scalar::Vy,
             layout.allocSize(HybridQuantity::Scalar::Vy)};
    Field Vz{"ions_bulkVel_z", HybridQuantity::Scalar::Vz,
             layout.allocSize(HybridQuantity::Scalar::Vz)};

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


    IonMoverTest()
        : ncells{100}
        , layout{{0.1}, {100u}, {{0.}}}
    {
        // First we need to set all buffer pointers of Electromag and Ions
        // Resources. Normally the ResourcesManger does that but here we don't
        // have a hierarchy etc. so we do that manually.


        EM.B.setBuffer("EM_B_x", &Bx);
        EM.B.setBuffer("EM_B_y", &By);
        EM.B.setBuffer("EM_B_z", &Bz);

        EM.E.setBuffer("EM_E_x", &Ex);
        EM.E.setBuffer("EM_E_y", &Ey);
        EM.E.setBuffer("EM_E_z", &Ez);

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

        protonPack.domainParticles        = &protonDomain;
        protonPack.patchGhostParticles    = &protonPatchGhost;
        protonPack.levelGhostParticles    = &protonLevelGhost;
        protonPack.levelGhostParticlesOld = &protonLevelGhostOld;
        protonPack.levelGhostParticlesNew = &protonLevelGhostNew;

        populations[0].setBuffer("ions_protons", &protonPack);

        populations[1].setBuffer("ions_alpha_rho", &alphaDensity);
        populations[1].flux().setBuffer("ions_alpha_flux_x", &alphaFx);
        populations[1].flux().setBuffer("ions_alpha_flux_y", &alphaFy);
        populations[1].flux().setBuffer("ions_alpha_flux_z", &alphaFz);

        alphaPack.domainParticles        = &alphaDomain;
        alphaPack.patchGhostParticles    = &alphaPatchGhost;
        alphaPack.levelGhostParticles    = &alphaLevelGhost;
        alphaPack.levelGhostParticlesOld = &alphaLevelGhostOld;
        alphaPack.levelGhostParticlesNew = &alphaLevelGhostNew;

        populations[1].setBuffer("ions_alpha", &alphaPack);


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
                // the number of ghost cells depends on the interpolator order
                if constexpr (interp_order == 1)
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
                                     return particle.iCell[0] == firstAMRCell[0];
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
                                     return particle.iCell[0] == lastAMRCell[0];
                                 });


                } // end first order
            }     // end 1D
        }         // end pop loop
    }
};



using DimInterps = ::testing::Types<DimInterp<1, 1>, DimInterp<1, 2>, DimInterp<1, 3>>;


TYPED_TEST_SUITE(IonMoverTest, DimInterps);


TYPED_TEST(IonMoverTest, bablabla)
{
    std::cout << 2 << "\n";
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
