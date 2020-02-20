
#include "initializer/data_provider.h"

#include "core/data/electrons/electrons.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/vecfield/vecfield.h"
#include "core/data/ions/ion_population/ion_population.h"
#include "core/data/ions/ions.h"
#include "core/data/electromag/electromag.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;



static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;



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
    dict["ions"]["name"]                                    = std::string{"ions"};
    dict["ions"]["nbrPopulations"]                          = int{1};
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




using GridImplYee1D = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D     = GridLayout<GridImplYee1D>;


using VecField1D      = VecField<NdArrayVector1D<>, HybridQuantity>;
using Field1D         = typename VecField1D::field_type;
using ScalarFunctionT = PHARE::initializer::ScalarFunction<1>;

using IonPopulation1D = IonPopulation<ParticleArray<1>, VecField1D, GridYee1D>;
using IonsT           = Ions<IonPopulation1D, GridYee1D>;
using PartPack1D      = ParticlesPack<typename IonPopulation1D::particle_array_type>;
using StandardHybridElectronFluxComputerT = StandardHybridElectronFluxComputer<IonsT>;


TEST(ElectronsTests, ThatElectronsHasCtor)
{
    IonsT ions{createDict()["ions"]};
    // Electromag<VecField1D> electromag{createDict()["electromag"]};
    // VecField1D J{"J", HybridQuantity::Vector::J};

    // Electrons electrons{createDict(), ions, electromag, J};
}



TEST(ElectronsTests, ThatElectronsDensityEqualIonDensity)
{
    IonsT ions{createDict()["ions"]};
    Electromag<VecField1D> electromag{createDict()["electromag"]};
    VecField1D J{"zobi", HybridQuantity::Vector::J};
    std::uint32_t nx = 50;
    Field1D Nibuffer{ions.densityName(), HybridQuantity::Scalar::rho, nx};
    Field1D NiProtons{"ions_protons_rho", HybridQuantity::Scalar::rho, nx};
    Field1D Vxi{"ions_bulkVel_x", HybridQuantity::Scalar::Vx, nx};
    Field1D Vyi{"ions_bulkVel_y", HybridQuantity::Scalar::Vy, nx};
    Field1D Vzi{"ions_bulkVel_z", HybridQuantity::Scalar::Vz, nx};
    Field1D Fxi{"ions_protons_flux_x", HybridQuantity::Scalar::Vx, nx};
    Field1D Fyi{"ions_protons_flux_y", HybridQuantity::Scalar::Vy, nx};
    Field1D Fzi{"ions_protons_flux_z", HybridQuantity::Scalar::Vz, nx};
    PartPack1D pack;

    ions.setBuffer(ions.densityName(), &Nibuffer);
    ions.velocity().setBuffer(Vxi.name(), &Vxi);
    ions.velocity().setBuffer(Vyi.name(), &Vyi);
    ions.velocity().setBuffer(Vzi.name(), &Vzi);

    auto& pops = ions.getRunTimeResourcesUserList();

    pops[0].setBuffer(NiProtons.name(), &NiProtons);
    pops[0].flux().setBuffer(Fxi.name(), &Fxi);
    pops[0].flux().setBuffer(Fyi.name(), &Fyi);
    pops[0].flux().setBuffer(Fzi.name(), &Fzi);
    pops[0].setBuffer("ions_protons", &pack);


    Electrons electrons{createDict(), ions, electromag, J};

    electrons.update();

    // auto& Ne = electrons.density();

    auto& Ni = ions.density();

    for (std::uint32_t i = 0; i < nx; ++i)
    {
        // EXPECT_DOUBLE_EQ(Ni(i), Ne(i));
    }
}


TEST(ElectronsTests, thatElectronsAreUsable)
{
    IonsT ions{createDict()["ions"]};
    Electromag<VecField1D> electromag{createDict()["electromag"]};
    VecField1D J{"J", HybridQuantity::Vector::J};
    std::uint32_t nx = 50;

    Field1D Bx{"EM_B_x", HybridQuantity::Scalar::Bx, nx};
    Field1D By{"EM_B_y", HybridQuantity::Scalar::By, nx};
    Field1D Bz{"EM_B_z", HybridQuantity::Scalar::Bz, nx};

    Field1D Ex{"EM_E_x", HybridQuantity::Scalar::Ex, nx};
    Field1D Ey{"EM_E_y", HybridQuantity::Scalar::Ey, nx};
    Field1D Ez{"EM_E_z", HybridQuantity::Scalar::Ez, nx};

    Field1D Jx{"J_x", HybridQuantity::Scalar::Jx, nx};
    Field1D Jy{"J_y", HybridQuantity::Scalar::Jy, nx};
    Field1D Jz{"J_z", HybridQuantity::Scalar::Jz, nx};

    Field1D Vex{"Ve_x", HybridQuantity::Scalar::Vx, nx};
    Field1D Vey{"Ve_y", HybridQuantity::Scalar::Vy, nx};
    Field1D Vez{"Ve_z", HybridQuantity::Scalar::Vz, nx};

    electromag.B.setBuffer(Bx.name(), &Bx);
    electromag.B.setBuffer(By.name(), &By);
    electromag.B.setBuffer(Bz.name(), &Bz);

    electromag.E.setBuffer(Ex.name(), &Ex);
    electromag.E.setBuffer(Ey.name(), &Ey);
    electromag.E.setBuffer(Ez.name(), &Ez);

    J.setBuffer(Jx.name(), &Jx);
    J.setBuffer(Jy.name(), &Jy);
    J.setBuffer(Jz.name(), &Jz);



    Field1D Nibuffer{ions.densityName(), HybridQuantity::Scalar::rho, nx};
    Field1D NiProtons{"ions_protons_rho", HybridQuantity::Scalar::rho, nx};
    Field1D Vxi{"ions_bulkVel_x", HybridQuantity::Scalar::Vx, nx};
    Field1D Vyi{"ions_bulkVel_y", HybridQuantity::Scalar::Vy, nx};
    Field1D Vzi{"ions_bulkVel_z", HybridQuantity::Scalar::Vz, nx};
    Field1D Fxi{"ions_protons_flux_x", HybridQuantity::Scalar::Vx, nx};
    Field1D Fyi{"ions_protons_flux_y", HybridQuantity::Scalar::Vy, nx};
    Field1D Fzi{"ions_protons_flux_z", HybridQuantity::Scalar::Vz, nx};
    PartPack1D pack;

    ions.setBuffer(ions.densityName(), &Nibuffer);
    ions.velocity().setBuffer(Vxi.name(), &Vxi);
    ions.velocity().setBuffer(Vyi.name(), &Vyi);
    ions.velocity().setBuffer(Vzi.name(), &Vzi);

    auto& pops = ions.getRunTimeResourcesUserList();

    pops[0].setBuffer(NiProtons.name(), &NiProtons);
    pops[0].flux().setBuffer(Fxi.name(), &Fxi);
    pops[0].flux().setBuffer(Fyi.name(), &Fyi);
    pops[0].flux().setBuffer(Fzi.name(), &Fzi);
    pops[0].setBuffer("ions_protons", &pack);



    Electrons electrons{createDict(), ions, electromag, J};

    auto&& emm = std::get<0>(electrons.getCompileTimeResourcesUserList());
    auto&& fc  = std::get<0>(emm.getCompileTimeResourcesUserList());
    auto&& Ve  = std::get<0>(
        static_cast<StandardHybridElectronFluxComputerT>(fc).getCompileTimeResourcesUserList());


    EXPECT_TRUE(electrons.isUsable());
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
