#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>

#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/ions/ion_population/particle_pack.h"

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

  - check the moments are equal to prescribed values
  - save moments and particles
  - apply tests

 */



double density(double x)
{
    return /*x * +*/ 2.;
}

double vx(double /*x*/)
{
    return 1.;
}


double vy(double /*x*/)
{
    return 1.;
}


double vz(double /*x*/)
{
    return 1.;
}


double vthx(double /*x*/)
{
    return 1.;
}


double vthy(double /*x*/)
{
    return 1.;
}


double vthz(double /*x*/)
{
    return 1.;
}


double bx(double x)
{
    return x /* + 1.*/;
}

double by(double x)
{
    return x /* + 2.*/;
}

double bz(double x)
{
    return x /*+ 3.*/;
}

double ex(double x)
{
    return x /* + 4.*/;
}

double ey(double x)
{
    return x /* + 5.*/;
}

double ez(double x)
{
    return x /* + 6.*/;
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

    ParticlesPack<ParticleArray> protonPack;
    ParticlesPack<ParticleArray> alphaPack;


    IonMoverTest()
        : ncells{100}
        , layout{{0.1}, {100u}, {{0.}}}
    {
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
        populations[0].setBuffer("ions_protons", &protonPack);

        populations[1].setBuffer("ions_alpha_rho", &alphaDensity);
        populations[1].flux().setBuffer("ions_alpha_flux_x", &alphaFx);
        populations[1].flux().setBuffer("ions_alpha_flux_y", &alphaFy);
        populations[1].flux().setBuffer("ions_alpha_flux_z", &alphaFz);
        populations[1].setBuffer("ions_alpha", &alphaPack);

        EM.initialize(layout);
    }
};



using DimInterps = ::testing::Types<DimInterp<1, 1>, DimInterp<1, 2>, DimInterp<1, 3>>;


TYPED_TEST_SUITE(IonMoverTest, DimInterps);


TYPED_TEST(IonMoverTest, bablabla)
{
    std::cout << 2 << "\n";
}


TEST(IonMoverTT, balbalbalbalablabla) {}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
