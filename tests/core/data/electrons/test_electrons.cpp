
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




class Electrons1DTest : public ::testing::Test
{
protected:
    IonsT ions;
    Electromag<VecField1D> electromag;
    VecField1D J;
    std::uint32_t nx = 50;
    Field1D Nibuffer;
    Field1D NiProtons;
    Field1D Vix;
    Field1D Viy;
    Field1D Viz;
    Field1D Fxi;
    Field1D Fyi;
    Field1D Fzi;
    PartPack1D pack;
    Field1D Vex;
    Field1D Vey;
    Field1D Vez;
    Field1D Jx;
    Field1D Jy;
    Field1D Jz;
    Field1D Bx;
    Field1D By;
    Field1D Bz;
    Field1D Ex;
    Field1D Ey;
    Field1D Ez;
    Electrons<IonsT, Electromag<VecField1D>> electrons;

public:
    Electrons1DTest()
        : ions{createDict()["ions"]}
        , electromag{createDict()["electromag"]}
        , J{"J", HybridQuantity::Vector::J}
        , Nibuffer{ions.densityName(), HybridQuantity::Scalar::rho, nx}
        , NiProtons{"ions_protons_rho", HybridQuantity::Scalar::rho, nx}
        , Vix{"ions_bulkVel_x", HybridQuantity::Scalar::Vx, nx}
        , Viy{"ions_bulkVel_y", HybridQuantity::Scalar::Vy, nx}
        , Viz{"ions_bulkVel_z", HybridQuantity::Scalar::Vz, nx}
        , Fxi{"ions_protons_flux_x", HybridQuantity::Scalar::Vx, nx}
        , Fyi{"ions_protons_flux_y", HybridQuantity::Scalar::Vy, nx}
        , Fzi{"ions_protons_flux_z", HybridQuantity::Scalar::Vz, nx}
        , Vex{"StandardHybridElectronFluxComputer_Ve_x", HybridQuantity::Scalar::Vx, nx}
        , Vey{"StandardHybridElectronFluxComputer_Ve_y", HybridQuantity::Scalar::Vy, nx}
        , Vez{"StandardHybridElectronFluxComputer_Ve_z", HybridQuantity::Scalar::Vz, nx}
        , Jx{"J_x", HybridQuantity::Scalar::Jx, nx}
        , Jy{"J_y", HybridQuantity::Scalar::Jy, nx}
        , Jz{"J_z", HybridQuantity::Scalar::Jz, nx}
        , Bx{"EM_B_x", HybridQuantity::Scalar::Bx, nx}
        , By{"EM_B_y", HybridQuantity::Scalar::By, nx}
        , Bz{"EM_B_z", HybridQuantity::Scalar::Bz, nx}
        , Ex{"EM_E_x", HybridQuantity::Scalar::Ex, nx}
        , Ey{"EM_E_y", HybridQuantity::Scalar::Ey, nx}
        , Ez{"EM_E_z", HybridQuantity::Scalar::Ez, nx}
        , electrons{createDict(), ions, electromag, J}
    {
        electromag.B.setBuffer(Bx.name(), &Bx);
        electromag.B.setBuffer(By.name(), &By);
        electromag.B.setBuffer(Bz.name(), &Bz);

        electromag.E.setBuffer(Ex.name(), &Ex);
        electromag.E.setBuffer(Ey.name(), &Ey);
        electromag.E.setBuffer(Ez.name(), &Ez);

        J.setBuffer(Jx.name(), &Jx);
        J.setBuffer(Jy.name(), &Jy);
        J.setBuffer(Jz.name(), &Jz);

        ions.setBuffer(ions.densityName(), &Nibuffer);
        ions.velocity().setBuffer(Vix.name(), &Vix);
        ions.velocity().setBuffer(Viy.name(), &Viy);
        ions.velocity().setBuffer(Viz.name(), &Viz);

        auto& pops = ions.getRunTimeResourcesUserList();

        pops[0].setBuffer(NiProtons.name(), &NiProtons);
        pops[0].flux().setBuffer(Fxi.name(), &Fxi);
        pops[0].flux().setBuffer(Fyi.name(), &Fyi);
        pops[0].flux().setBuffer(Fzi.name(), &Fzi);
        pops[0].setBuffer("ions_protons", &pack);

        auto&& emm = std::get<0>(electrons.getCompileTimeResourcesUserList());
        auto&& fc  = std::get<0>(emm.getCompileTimeResourcesUserList());
        auto&& Ve  = std::get<0>(fc.getCompileTimeResourcesUserList());

        Ve.setBuffer(Vex.name(), &Vex);
        Ve.setBuffer(Vey.name(), &Vey);
        Ve.setBuffer(Vez.name(), &Vez);
    }
};



TEST_F(Electrons1DTest, ThatElectronsHasCtor) {}



TEST_F(Electrons1DTest, ThatElectronsAreUsable)
{
    EXPECT_TRUE(electrons.isUsable());
}




TEST_F(Electrons1DTest, ThatElectronsDensityEqualIonDensity)
{
    electrons.update();

    auto& Ne = electrons.density();
    auto& Ni = ions.density();

    for (std::uint32_t i = 0; i < nx; ++i)
    {
        EXPECT_DOUBLE_EQ(Ni(i), Ne(i));
    }
}


TEST_F(Electrons1DTest, ThatElectronsVelocityEqualIonVelocityMinusJ)
{
    electrons.update();

    auto& Ne = electrons.density();

    for (std::uint32_t i = 0; i < nx; ++i)
    {
        // EXPECT_DOUBLE_EQ(Vex(i), Vix(i) - Jx(i) / Ne(i));
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
