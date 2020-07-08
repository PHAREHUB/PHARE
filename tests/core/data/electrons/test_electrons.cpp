
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

const double Te = 0.12;


using ScalarFunctionT = PHARE::initializer::ScalarFunction<1>;

PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;
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

    dict["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    dict["electrons"]["pressure_closure"]["Te"]   = Te;

    return dict;
}



using GridImplYee = GridLayoutImplYee<dim, interpOrder>;
using GridYee     = GridLayout<GridImplYee>;

using VecFieldND      = VecField<NdArrayVector<dim>, HybridQuantity>;
using FieldND         = typename VecFieldND::field_type;
using ScalarFunctionT = PHARE::initializer::ScalarFunction<dim>;

using IonPopulationND = IonPopulation<ParticleArray<dim>, VecFieldND, GridYee>;
using IonsT           = Ions<IonPopulationND, GridYee>;
using PartPackND      = ParticlesPack<typename IonPopulationND::particle_array_type>;
using StandardHybridElectronFluxComputerT = StandardHybridElectronFluxComputer<IonsT>;




// https://stackoverflow.com/questions/46101569/compile-time-constructor-switch-in-c
template<int>
struct theDim
{
};


class nDLayout
{
    constexpr nDLayout(std::array<double, dim> mesh, std::array<uint32, dim> numofcells,
                       Point<double, dim> origin)
        : meshSize{mesh}
        , nbrCells{numofcells}
        , origin{origin}
    {
    }

    constexpr nDLayout(theDim<1>)
        : nDLayout{{{0.1}}, {{50}}, Point<double, 1>{0.}}
    {
    }

    // constexpr nDLayout(theDim<2>)
    //    : nDLayout{{{0.1, 0.2}}, {{50, 30}}, {0., 0.}}
    //{
    //}

private:
    std::array<double, dim> meshSize;
    std::array<uint32, dim> nbrCells;
    Point<double, dim> origin;
};



class ElectronsTest : public ::testing::Test
{
protected:
    GridYee layout;
    IonsT ions;
    Electromag<VecFieldND> electromag;
    VecFieldND J;
    FieldND Nibuffer;
    FieldND NiProtons;
    FieldND Vix;
    FieldND Viy;
    FieldND Viz;
    FieldND Fxi;
    FieldND Fyi;
    FieldND Fzi;
    PartPackND pack;
    FieldND Vex;
    FieldND Vey;
    FieldND Vez;
    FieldND Jx;
    FieldND Jy;
    FieldND Jz;
    Electrons<IonsT> electrons;
    FieldND Pe;

public:
    ElectronsTest()
        : // layout{nDLayout(theDim<dim>)}
        layout{{{0.1}}, {{50}}, Point<double, dim>{0.}}
        , ions{createDict()["ions"]}
        , electromag{createDict()["electromag"]}
        , J{"J", HybridQuantity::Vector::J}
        , Nibuffer{ions.densityName(), HybridQuantity::Scalar::rho,
                   layout.allocSize(HybridQuantity::Scalar::rho)}
        , NiProtons{"protons_rho", HybridQuantity::Scalar::rho,
                    layout.allocSize(HybridQuantity::Scalar::rho)}
        , Vix{"bulkVel_x", HybridQuantity::Scalar::Vx, layout.allocSize(HybridQuantity::Scalar::Vx)}
        , Viy{"bulkVel_y", HybridQuantity::Scalar::Vy, layout.allocSize(HybridQuantity::Scalar::Vy)}
        , Viz{"bulkVel_z", HybridQuantity::Scalar::Vz, layout.allocSize(HybridQuantity::Scalar::Vz)}
        , Fxi{"protons_flux_x", HybridQuantity::Scalar::Vx,
              layout.allocSize(HybridQuantity::Scalar::Vx)}
        , Fyi{"protons_flux_y", HybridQuantity::Scalar::Vy,
              layout.allocSize(HybridQuantity::Scalar::Vy)}
        , Fzi{"protons_flux_z", HybridQuantity::Scalar::Vz,
              layout.allocSize(HybridQuantity::Scalar::Vz)}
        , Vex{"StandardHybridElectronFluxComputer_Ve_x", HybridQuantity::Scalar::Vx,
              layout.allocSize(HybridQuantity::Scalar::Vx)}
        , Vey{"StandardHybridElectronFluxComputer_Ve_y", HybridQuantity::Scalar::Vy,
              layout.allocSize(HybridQuantity::Scalar::Vy)}
        , Vez{"StandardHybridElectronFluxComputer_Ve_z", HybridQuantity::Scalar::Vz,
              layout.allocSize(HybridQuantity::Scalar::Vz)}
        , Jx{"J_x", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"J_y", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"J_z", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , electrons{createDict()["electrons"], ions, J}
        , Pe{"Pe", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)}
    {
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
        pops[0].setBuffer("protons", &pack);

        auto&& emm = std::get<0>(electrons.getCompileTimeResourcesUserList());
        auto&& fc  = std::get<0>(emm.getCompileTimeResourcesUserList());
        auto&& Ve  = std::get<0>(fc.getCompileTimeResourcesUserList());

        Ve.setBuffer(Vex.name(), &Vex);
        Ve.setBuffer(Vey.name(), &Vey);
        Ve.setBuffer(Vez.name(), &Vez);

        auto&& pc = std::get<1>(emm.getCompileTimeResourcesUserList());

        pc.setBuffer(Pe.name(), &Pe);

        if constexpr (dim == 1)
        {
            auto fill = [this](FieldND& field, auto const& filler) {
                auto gsi_X = this->layout.ghostStartIndex(field, Direction::X);
                auto gei_X = this->layout.ghostEndIndex(field, Direction::X);

                for (auto ix = gsi_X; ix <= gei_X; ++ix)
                {
                    auto point = this->layout.fieldNodeCoordinates(field, Point<double, 1>{0.}, ix);
                    field(ix)  = filler(point[0]);
                }
            };

            fill(Vix, [](double x) { return std::cosh(0.2 * x); });
            fill(Viy, [](double x) { return std::cosh(0.3 * x); });
            fill(Viz, [](double x) { return std::cosh(0.4 * x); });

            fill(Jx, [](double x) { return std::sinh(0.2 * x); });
            fill(Jy, [](double x) { return std::sinh(0.3 * x); });
            fill(Jz, [](double x) { return std::sinh(0.4 * x); });

            fill(Nibuffer, [](double x) { return std::cosh(0.1 * x); });
        }
        else if constexpr (dim == 2)
        {
        }
        else if constexpr (dim == 3)
        {
        }
    }


    void TearDown() override
    {
        J.setBuffer(Jx.name(), nullptr);
        J.setBuffer(Jy.name(), nullptr);
        J.setBuffer(Jz.name(), nullptr);

        ions.setBuffer(ions.densityName(), nullptr);
        ions.velocity().setBuffer(Vix.name(), nullptr);
        ions.velocity().setBuffer(Viy.name(), nullptr);
        ions.velocity().setBuffer(Viz.name(), nullptr);

        auto& pops = ions.getRunTimeResourcesUserList();

        pops[0].setBuffer(NiProtons.name(), static_cast<FieldND*>(nullptr));
        pops[0].flux().setBuffer(Fxi.name(), nullptr);
        pops[0].flux().setBuffer(Fyi.name(), nullptr);
        pops[0].flux().setBuffer(Fzi.name(), nullptr);
        pops[0].setBuffer("protons", static_cast<PartPackND*>(nullptr));

        auto&& emm = std::get<0>(electrons.getCompileTimeResourcesUserList());
        auto&& fc  = std::get<0>(emm.getCompileTimeResourcesUserList());
        auto&& Ve  = std::get<0>(fc.getCompileTimeResourcesUserList());

        Ve.setBuffer(Vex.name(), nullptr);
        Ve.setBuffer(Vey.name(), nullptr);
        Ve.setBuffer(Vez.name(), nullptr);

        auto&& pc = std::get<1>(emm.getCompileTimeResourcesUserList());

        pc.setBuffer(Pe.name(), nullptr);
    }
};



TEST_F(ElectronsTest, ThatElectronsHasCtor) {}



TEST_F(ElectronsTest, ThatElectronsAreUsable)
{
    EXPECT_TRUE(electrons.isUsable());
}




TEST_F(ElectronsTest, ThatElectronsDensityEqualIonDensity)
{
    electrons.update(layout);

    auto& Ne = electrons.density();
    auto& Ni = ions.density();

    if constexpr (dim == 1)
    {
        auto psi_X = this->layout.physicalStartIndex(Ne, Direction::X);
        auto pei_X = this->layout.physicalEndIndex(Ne, Direction::X);

        for (std::uint32_t i = psi_X; i < pei_X; ++i)
        {
            EXPECT_DOUBLE_EQ(Ni(i), Ne(i));
        }
    }
    else if constexpr (dim == 2)
    {
    }
    else if constexpr (dim == 3)
    {
    }
}


TEST_F(ElectronsTest, ThatElectronsVelocityEqualIonVelocityMinusJ)
{
    electrons.update(layout);

    auto& Ne = electrons.density();


    auto check = [this](FieldND const& Vecomp, FieldND const& Vicomp, FieldND const& Jcomp,
                        FieldND const& Ne_, auto const& projector) {
        if constexpr (dim == 1)
        {
            auto psi_X = this->layout.physicalStartIndex(Vicomp, Direction::X);
            auto pei_X = this->layout.physicalEndIndex(Vicomp, Direction::X);

            for (std::uint32_t i = psi_X; i < pei_X; ++i)
            {
                auto const JOnV = GridYee::project(Jcomp, {i}, projector());

                EXPECT_DOUBLE_EQ(Vecomp(i), Vicomp(i) - JOnV / Ne_(i));
            }
        }
        else if constexpr (dim == 2)
        {
        }
        else if constexpr (dim == 3)
        {
        }
    };

    check(Vex, Vix, Jx, Ne, GridYee::JxToMoments);
    check(Vey, Viy, Jy, Ne, GridYee::JyToMoments);
    check(Vez, Viz, Jz, Ne, GridYee::JzToMoments);
}



TEST_F(ElectronsTest, ThatElectronsPressureEqualsNeTe)
{
    electrons.update(layout);

    auto& Ne_ = electrons.density();
    auto& Pe_ = electrons.pressure();

    if constexpr (dim == 1)
    {
        auto psi_X = this->layout.physicalStartIndex(Ne_, Direction::X);
        auto pei_X = this->layout.physicalEndIndex(Ne_, Direction::X);

        for (std::uint32_t i = psi_X; i < pei_X; ++i)
        {
            EXPECT_DOUBLE_EQ(Pe_(i), Ne_(i) * Te);
        }
    }
    else if constexpr (dim == 2)
    {
    }
    else if constexpr (dim == 3)
    {
    }
}



template<std::size_t dim, std::size_t interpO>
struct dimAndInterpOrder
{
    static constexpr auto dimension = dim;
    static constexpr auto interp    = interpO;
};

using MyDimAndInterpOrders
    = ::testing::Types<dimAndInterpOrder<1, 1>, dimAndInterpOrder<1, 2>, dimAndInterpOrder<1, 3>>;

// typedef ::testing::Types<Dimension<1>,Dimension<2>,Dimension<3>




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
