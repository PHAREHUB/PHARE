

#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "initializer/data_provider.hpp"

#include "core/data/grid/grid.hpp"
#include "core/data/ions/ions.hpp"
#include "src/core/utilities/types.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/electrons/electrons.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"

#include "gtest/gtest.h"

#include "tests/initializer/init_functions.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"


#include <string_view>

using namespace PHARE::core;

double const Te = 0.12;



template<int dim>
void setPopulation(PHARE::initializer::PHAREDict& dict, std::string const& popKey,
                   std::string const& name, double const mass, double const charge)
{
    using InitFunctionT = PHARE::initializer::InitFunction<dim>;

    auto density = makeSharedPtr<dim>();
    auto vx      = makeSharedPtr<dim>();
    auto vy      = makeSharedPtr<dim>();
    auto vz      = makeSharedPtr<dim>();
    auto vthx    = makeSharedPtr<dim>();
    auto vthy    = makeSharedPtr<dim>();
    auto vthz    = makeSharedPtr<dim>();

    auto& pop                           = dict["ions"][popKey];
    pop["name"]                         = name;
    pop["mass"]                         = mass;
    pop["particle_initializer"]["name"] = std::string{"maxwellian"};

    pop["particle_initializer"]["density"]            = static_cast<InitFunctionT>(density);
    pop["particle_initializer"]["bulk_velocity_x"]    = static_cast<InitFunctionT>(vx);
    pop["particle_initializer"]["bulk_velocity_y"]    = static_cast<InitFunctionT>(vy);
    pop["particle_initializer"]["bulk_velocity_z"]    = static_cast<InitFunctionT>(vz);
    pop["particle_initializer"]["thermal_velocity_x"] = static_cast<InitFunctionT>(vthx);
    pop["particle_initializer"]["thermal_velocity_y"] = static_cast<InitFunctionT>(vthy);
    pop["particle_initializer"]["thermal_velocity_z"] = static_cast<InitFunctionT>(vthz);

    pop["particle_initializer"]["nbr_part_per_cell"] = int{100};
    pop["particle_initializer"]["charge"]            = charge;
    pop["particle_initializer"]["basis"]             = std::string{"cartesian"};
}


template<int dim>
PHARE::initializer::PHAREDict createDict()
{
    using InitFunctionT = PHARE::initializer::InitFunction<dim>;

    auto bx = makeSharedPtr<dim>();
    auto by = makeSharedPtr<dim>();
    auto bz = makeSharedPtr<dim>();
    auto ex = makeSharedPtr<dim>();
    auto ey = makeSharedPtr<dim>();
    auto ez = makeSharedPtr<dim>();


    PHARE::initializer::PHAREDict dict;
    dict["ions"]["nbrPopulations"] = std::size_t{1};
    setPopulation<dim>(dict, "pop0", "protons", /*mass=*/1., /*charge=*/-1.);

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["electric"]["initializer"]["x_component"] = static_cast<InitFunctionT>(ex);
    dict["electromag"]["electric"]["initializer"]["y_component"] = static_cast<InitFunctionT>(ey);
    dict["electromag"]["electric"]["initializer"]["z_component"] = static_cast<InitFunctionT>(ez);

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<InitFunctionT>(bx);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<InitFunctionT>(by);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<InitFunctionT>(bz);

    dict["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    dict["electrons"]["pressure_closure"]["Te"]   = Te;

    return dict;
}


template<int dim, int interp>
class NDlayout
{
    NDlayout() {}

    using nDL = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;

public:
    static nDL create()
    {
        if constexpr (dim == 1)
        {
            return {{{0.1}}, {{50}}, {0.}};
        }
        else if constexpr (dim == 2)
        {
            return {{{0.1, 0.2}}, {{50, 40}}, {0., 0.}};
        }
        else if constexpr (dim == 3)
        {
            return {{{0.1, 0.2, 0.3}}, {{50, 40, 30}}, {0., 0., 0.}};
        }
    }
};



template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct ElectronsTest : public ::testing::Test
{
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();


    using GridYee = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;

    using GridND           = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using FieldND          = Field<dim, HybridQuantity::Scalar>;
    using VecFieldND       = VecField<FieldND, HybridQuantity>;
    using SymTensorFieldND = SymTensorField<FieldND, HybridQuantity>;
    using ParticleArray_t  = ParticleArray<dim>;
    using IonPopulationND  = IonPopulation<ParticleArray_t, VecFieldND, SymTensorFieldND>;
    using IonsT            = Ions<IonPopulationND, GridYee>;
    using PartPackND       = ParticlesPack<ParticleArray_t>;
    using StandardHybridElectronFluxComputerT = StandardHybridElectronFluxComputer<IonsT>;


    GridYee layout = NDlayout<dim, interp>::create();

    Electromag<VecFieldND> electromag;

    UsableVecField<dim> J, F, Ve, Vi, chargeFlux;
    UsableTensorField<dim> ionTensor, protonTensor;

    GridND ionChargeDensity, ionMassDensity, protonParticleDensity, protonChargeDensity, Pe;

    ParticleArray_t domainParticles{layout.AMRBox()};
    ParticleArray_t patchGhostParticles = domainParticles;
    ParticleArray_t levelGhostParticles = domainParticles;
    PartPackND pack{"particles", &domainParticles, &patchGhostParticles, &levelGhostParticles};

    IonsT ions;
    Electrons<IonsT> electrons;

    template<typename... Args>
    auto static _ions(Args&... args)
    {
        auto const& [ionFlux, ionChargeDensity, ionMassDensity, protonParticleDensity,
                     protonChargeDensity, Vi, chargeFlux, ionTensor, protonTensor, pack]
            = std::forward_as_tuple(args...);
        IonsT ions{createDict<dim>()["ions"]};
        {
            auto const& [V, m, d_c, d_m, c_f] = ions.getCompileTimeResourcesViewList();
            d_c.setBuffer(&ionChargeDensity);
            d_m.setBuffer(&ionMassDensity);
            Vi.set_on(V);
            chargeFlux.set_on(c_f);
            ionTensor.set_on(m);
        }
        auto& pops = ions.getRunTimeResourcesViewList();
        assert(pops.size() == 1);

        auto const& [F, m, Np, Nc, poppack] = pops[0].getCompileTimeResourcesViewList();
        Np.setBuffer(&protonParticleDensity);
        Nc.setBuffer(&protonChargeDensity);
        ionFlux.set_on(F);
        protonTensor.set_on(m);
        poppack.setBuffer(&pack);
        return ions;
    }


    ElectronsTest()
        : electromag{createDict<dim>()["electromag"]}
        , J{"J", layout, HybridQuantity::Vector::J}
        , F{"protons_flux", layout, HybridQuantity::Vector::V}
        , Ve{"StandardHybridElectronFluxComputer_Ve", layout, HybridQuantity::Vector::V}
        , Vi{"bulkVel", layout, HybridQuantity::Vector::V}
        , chargeFlux{"chargeFlux", layout, HybridQuantity::Vector::V}
        , ionTensor{"momentumTensor", layout, HybridQuantity::Tensor::M}
        , protonTensor{"protons_momentumTensor", layout, HybridQuantity::Tensor::M}
        , ionChargeDensity{"chargeDensity", HybridQuantity::Scalar::rho,
                           layout.allocSize(HybridQuantity::Scalar::rho)}
        , ionMassDensity{"massDensity", HybridQuantity::Scalar::rho,
                         layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonParticleDensity{"protons_particleDensity", HybridQuantity::Scalar::rho,
                                layout.allocSize(HybridQuantity::Scalar::rho)}
        , protonChargeDensity{"protons_chargeDensity", HybridQuantity::Scalar::rho,
                              layout.allocSize(HybridQuantity::Scalar::rho)}
        , Pe{"Pe", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)}
        , ions{_ions(F, ionChargeDensity, ionMassDensity, protonParticleDensity,
                     protonChargeDensity, Vi, chargeFlux, ionTensor, protonTensor, pack)}
        , electrons{createDict<dim>()["electrons"], ions, J}
    {
        auto&& emm = std::get<0>(electrons.getCompileTimeResourcesViewList());
        auto&& fc  = std::get<0>(emm.getCompileTimeResourcesViewList());


        Ve.set_on(std::get<0>(fc.getCompileTimeResourcesViewList()));


        auto&& pc          = std::get<1>(emm.getCompileTimeResourcesViewList());
        auto const& [_, P] = pc.getCompileTimeResourcesViewList();
        P.setBuffer(&Pe);

        auto const& [Jx, Jy, Jz] = J();
        auto const& [Fx, Fy, Fz] = chargeFlux();


        if constexpr (dim == 1)
        {
            auto fill = [this](auto& field, auto const& filler) {
                layout.ghost_amr_lcl_idx(field, [&](auto const& amr_idx, auto const& lcl_idx) {
                    auto const point = this->layout.fieldNodeCoordinates(field, amr_idx);
                    field(lcl_idx)   = filler(point[0]);
                });
            };

            fill(Fx, [](double x) { return std::cosh(0.2 * x); });
            fill(Fy, [](double x) { return std::cosh(0.3 * x); });
            fill(Fz, [](double x) { return std::cosh(0.4 * x); });

            fill(Jx, [](double x) { return std::sinh(0.2 * x); });
            fill(Jy, [](double x) { return std::sinh(0.3 * x); });
            fill(Jz, [](double x) { return std::sinh(0.4 * x); });

            fill(ionChargeDensity, [](double x) { return std::cosh(0.1 * x); });
        }
        else if constexpr (dim == 2)
        {
            auto fill = [this](auto& field, auto const& filler) {
                for_each_point(layout.AMRGhostBoxFor(field), [&](auto const& amr_idx) {
                    auto point = this->layout.fieldNodeCoordinates(field, amr_idx);
                    field(layout.AMRToLocal(amr_idx)) = filler(point[0], point[1]);
                });
            };

            fill(Fx, [](double x, double y) { return std::cosh(0.2 * x) * std::cosh(0.2 * y); });
            fill(Fy, [](double x, double y) { return std::cosh(0.3 * x) * std::cosh(0.3 * y); });
            fill(Fz, [](double x, double y) { return std::cosh(0.4 * x) * std::cosh(0.4 * y); });

            fill(Jx, [](double x, double y) { return std::sinh(0.2 * x) * std::sinh(0.2 * y); });
            fill(Jy, [](double x, double y) { return std::sinh(0.3 * x) * std::sinh(0.3 * y); });
            fill(Jz, [](double x, double y) { return std::sinh(0.4 * x) * std::sinh(0.4 * y); });

            fill(ionChargeDensity,
                 [](double x, double y) { return std::cosh(0.1 * x) * std::cosh(0.1 * y); });
        }
        else if constexpr (dim == 3)
        {
            auto fill = [this](auto& field, auto const& filler) {
                for_each_point(layout.AMRGhostBoxFor(field), [&](auto const& amr_idx) {
                    auto point = this->layout.fieldNodeCoordinates(field, amr_idx);
                    field(layout.AMRToLocal(amr_idx)) = filler(point[0], point[1], point[2]);
                });
            };

            fill(Fx, [](double x, double y, double z) {
                return std::cosh(0.2 * x) * std::cosh(0.2 * y) * std::cosh(0.2 * z);
            });
            fill(Fy, [](double x, double y, double z) {
                return std::cosh(0.3 * x) * std::cosh(0.3 * y) * std::cosh(0.3 * z);
            });
            fill(Fz, [](double x, double y, double z) {
                return std::cosh(0.4 * x) * std::cosh(0.4 * y) * std::cosh(0.4 * z);
            });

            fill(Jx, [](double x, double y, double z) {
                return std::sinh(0.2 * x) * std::sinh(0.2 * y) * std::sinh(0.2 * z);
            });
            fill(Jy, [](double x, double y, double z) {
                return std::sinh(0.3 * x) * std::sinh(0.3 * y) * std::sinh(0.3 * z);
            });
            fill(Jz, [](double x, double y, double z) {
                return std::sinh(0.4 * x) * std::sinh(0.4 * y) * std::sinh(0.4 * z);
            });

            fill(ionChargeDensity, [](double x, double y, double z) {
                return std::cosh(0.1 * x) * std::cosh(0.1 * y) * std::cosh(0.1 * z);
            });
        }
    }
};


using ElectronsTupleInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<1>, InterpConst<2>>,
                     std::pair<DimConst<1>, InterpConst<3>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<2>, InterpConst<2>>, std::pair<DimConst<2>, InterpConst<3>>,
                     std::pair<DimConst<3>, InterpConst<1>>, std::pair<DimConst<3>, InterpConst<2>>,
                     std::pair<DimConst<3>, InterpConst<3>>>;

TYPED_TEST_SUITE(ElectronsTest, ElectronsTupleInfos, );



TYPED_TEST(ElectronsTest, ThatElectronsHasCtor)
{
    static constexpr auto dim    = typename TypeParam::first_type{}();
    static constexpr auto interp = typename TypeParam::second_type{}();

    std::cout << __FILE__ << " " << __LINE__ << " " << dim << std::endl;
    std::cout << __FILE__ << " " << __LINE__ << " " << interp << std::endl;
}



TYPED_TEST(ElectronsTest, ThatElectronsAreUsable)
{
    auto& electrons = this->electrons;
    EXPECT_TRUE(electrons.isUsable());
}




TYPED_TEST(ElectronsTest, ThatElectronsDensityEqualIonDensity)
{
    static constexpr auto dim = typename TypeParam::first_type{}();

    auto& electrons = this->electrons;
    auto& layout    = this->layout;
    auto& ions      = this->ions;

    electrons.update(layout);

    auto& Ne = electrons.density();
    auto& Ni = ions.chargeDensity();

    if constexpr (dim == 1)
    {
        auto psi_X = layout.physicalStartIndex(Ne, Direction::X);
        auto pei_X = layout.physicalEndIndex(Ne, Direction::X);

        for (std::uint32_t i = psi_X; i < pei_X; ++i)
        {
            EXPECT_DOUBLE_EQ(Ni(i), Ne(i));
        }
    }
    else if constexpr (dim == 2)
    {
        auto psi_X = layout.physicalStartIndex(Ne, Direction::X);
        auto pei_X = layout.physicalEndIndex(Ne, Direction::X);
        auto psi_Y = layout.physicalStartIndex(Ne, Direction::Y);
        auto pei_Y = layout.physicalEndIndex(Ne, Direction::Y);

        for (std::uint32_t i = psi_X; i < pei_X; ++i)
        {
            for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
            {
                EXPECT_DOUBLE_EQ(Ni(i, j), Ne(i, j));
            }
        }
    }
    else if constexpr (dim == 3)
    {
        auto psi_X = layout.physicalStartIndex(Ne, Direction::X);
        auto pei_X = layout.physicalEndIndex(Ne, Direction::X);
        auto psi_Y = layout.physicalStartIndex(Ne, Direction::Y);
        auto pei_Y = layout.physicalEndIndex(Ne, Direction::Y);
        auto psi_Z = layout.physicalStartIndex(Ne, Direction::Z);
        auto pei_Z = layout.physicalEndIndex(Ne, Direction::Z);

        for (std::uint32_t i = psi_X; i < pei_X; ++i)
        {
            for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
            {
                for (std::uint32_t k = psi_Z; k < pei_Z; ++k)
                {
                    EXPECT_DOUBLE_EQ(Ni(i, j, k), Ne(i, j, k));
                }
            }
        }
    }
}


TYPED_TEST(ElectronsTest, ThatElectronsVelocityEqualChargeFluxMinusJOverNe)
{
    static constexpr auto dim    = typename TypeParam::first_type{}();
    static constexpr auto interp = typename TypeParam::second_type{}();
    using GridYee   = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;
    auto& electrons = this->electrons;
    auto& layout    = this->layout;
    electrons.update(layout);
    auto& Ne = electrons.density();

    auto check = [&layout]<auto projector>(auto const& Vecomp, auto const& Fcomp, auto const& Jcomp,
                                           auto const& Ne_) -> void {
        if constexpr (dim == 1)
        {
            auto psi_X = layout.physicalStartIndex(Fcomp, Direction::X);
            auto pei_X = layout.physicalEndIndex(Fcomp, Direction::X);
            for (std::uint32_t i = psi_X; i < pei_X; ++i)
            {
                auto const JOnV = GridYee::template project<projector>(Jcomp, {i});
                EXPECT_DOUBLE_EQ(Vecomp(i), (Fcomp(i) - JOnV) / Ne_(i));
            }
        }
        else if constexpr (dim == 2)
        {
            auto psi_X = layout.physicalStartIndex(Fcomp, Direction::X);
            auto pei_X = layout.physicalEndIndex(Fcomp, Direction::X);
            auto psi_Y = layout.physicalStartIndex(Fcomp, Direction::Y);
            auto pei_Y = layout.physicalEndIndex(Fcomp, Direction::Y);
            for (std::uint32_t i = psi_X; i < pei_X; ++i)
                for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
                {
                    auto const JOnV = GridYee::template project<projector>(Jcomp, {i, j});
                    EXPECT_DOUBLE_EQ(Vecomp(i, j), (Fcomp(i, j) - JOnV) / Ne_(i, j));
                }
        }
        else if constexpr (dim == 3)
        {
            auto psi_X = layout.physicalStartIndex(Fcomp, Direction::X);
            auto pei_X = layout.physicalEndIndex(Fcomp, Direction::X);
            auto psi_Y = layout.physicalStartIndex(Fcomp, Direction::Y);
            auto pei_Y = layout.physicalEndIndex(Fcomp, Direction::Y);
            auto psi_Z = layout.physicalStartIndex(Fcomp, Direction::Z);
            auto pei_Z = layout.physicalEndIndex(Fcomp, Direction::Z);
            for (std::uint32_t i = psi_X; i < pei_X; ++i)
                for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
                    for (std::uint32_t k = psi_Z; k < pei_Z; ++k)
                    {
                        auto const JOnV = GridYee::template project<projector>(Jcomp, {i, j, k});
                        EXPECT_DOUBLE_EQ(Vecomp(i, j, k), (Fcomp(i, j, k) - JOnV) / Ne_(i, j, k));
                    }
        }
    };

    auto const& [Jx, Jy, Jz]    = this->J();
    auto const& [Fx, Fy, Fz]    = this->chargeFlux();
    auto const& [Vex, Vey, Vez] = this->Ve();

    check.template operator()<GridYee::JxToMoments>(Vex, Fx, Jx, Ne);
    check.template operator()<GridYee::JyToMoments>(Vey, Fy, Jy, Ne);
    check.template operator()<GridYee::JzToMoments>(Vez, Fz, Jz, Ne);
}


TYPED_TEST(ElectronsTest, ThatElectronsPressureEqualsNeTe)
{
    static constexpr auto dim = typename TypeParam::first_type{}();

    auto& electrons = this->electrons;
    auto& layout    = this->layout;

    electrons.update(layout);

    auto& Ne_ = electrons.density();
    auto& Pe_ = electrons.pressure();

    if constexpr (dim == 1)
    {
        auto psi_X = layout.physicalStartIndex(Ne_, Direction::X);
        auto pei_X = layout.physicalEndIndex(Ne_, Direction::X);

        for (std::uint32_t i = psi_X; i < pei_X; ++i)
        {
            EXPECT_DOUBLE_EQ(Pe_(i), Ne_(i) * Te);
        }
    }
    else if constexpr (dim == 2)
    {
        auto psi_X = layout.physicalStartIndex(Ne_, Direction::X);
        auto pei_X = layout.physicalEndIndex(Ne_, Direction::X);
        auto psi_Y = layout.physicalStartIndex(Ne_, Direction::Y);
        auto pei_Y = layout.physicalEndIndex(Ne_, Direction::Y);

        for (std::uint32_t i = psi_X; i < pei_X; ++i)
        {
            for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
            {
                EXPECT_DOUBLE_EQ(Pe_(i, j), Ne_(i, j) * Te);
            }
        }
    }
    else if constexpr (dim == 3)
    {
        auto psi_X = layout.physicalStartIndex(Ne_, Direction::X);
        auto pei_X = layout.physicalEndIndex(Ne_, Direction::X);
        auto psi_Y = layout.physicalStartIndex(Ne_, Direction::Y);
        auto pei_Y = layout.physicalEndIndex(Ne_, Direction::Y);
        auto psi_Z = layout.physicalStartIndex(Ne_, Direction::Z);
        auto pei_Z = layout.physicalEndIndex(Ne_, Direction::Z);

        for (std::uint32_t i = psi_X; i < pei_X; ++i)
        {
            for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
            {
                for (std::uint32_t k = psi_Z; k < pei_Z; ++k)
                {
                    EXPECT_DOUBLE_EQ(Pe_(i, j, k), Ne_(i, j, k) * Te);
                }
            }
        }
    }
}




// two populations with different q/m drifting relative to each other:
// the electron velocity must be (sum_s q_s F_s - J) / Ne, see issue #1466
template<typename TypeInfo>
struct TwoPopulationElectronsTest : public ::testing::Test
{
    using Base = ElectronsTest<TypeInfo>;

    static constexpr auto dim = Base::dim;

    using GridYee         = typename Base::GridYee;
    using GridND          = typename Base::GridND;
    using ParticleArray_t = typename Base::ParticleArray_t;
    using PartPackND      = typename Base::PartPackND;
    using IonsT           = typename Base::IonsT;

    struct Population
    {
        std::string name;
        double mass, charge, density;
        std::array<double, 3> velocity;
    };

    static constexpr std::size_t nbrPops = 2;
    std::array<Population, nbrPops> const pops{{
        {"protons", /*mass=*/1., /*charge=*/1., /*density=*/1., {1., 0.5, -0.3}},
        {"sodium", /*mass=*/23., /*charge=*/2., /*density=*/0.1, {-0.2, 0.4, 0.1}},
    }};
    std::array<double, 3> const current{0.1, -0.2, 0.3};

    struct PopulationBuffers
    {
        PopulationBuffers(std::string const& name, GridYee const& layout)
            : flux{name + "_flux", layout, HybridQuantity::Vector::V}
            , momentumTensor{name + "_momentumTensor", layout, HybridQuantity::Tensor::M}
            , particleDensity{name + "_particleDensity", HybridQuantity::Scalar::rho,
                              layout.allocSize(HybridQuantity::Scalar::rho)}
            , chargeDensity{name + "_chargeDensity", HybridQuantity::Scalar::rho,
                            layout.allocSize(HybridQuantity::Scalar::rho)}
            , domainParticles{layout.AMRBox()}
            , pack{name, &domainParticles, &patchGhostParticles, &levelGhostParticles}
        {
        }

        UsableVecField<dim> flux;
        UsableTensorField<dim> momentumTensor;
        GridND particleDensity, chargeDensity;
        ParticleArray_t domainParticles;
        ParticleArray_t patchGhostParticles = domainParticles;
        ParticleArray_t levelGhostParticles = domainParticles;
        PartPackND pack;
    };


    static auto createTwoPopulationDict(auto const& pops_)
    {
        auto dict                      = createDict<dim>();
        dict["ions"]["nbrPopulations"] = nbrPops;
        for (std::size_t i = 0; i < nbrPops; ++i)
            setPopulation<dim>(dict, "pop" + std::to_string(i), pops_[i].name, pops_[i].mass,
                               pops_[i].charge);
        return dict;
    }


    IonsT makeIons()
    {
        IonsT ions_{dict["ions"]};
        {
            auto const& [V, m, d_c, d_m, c_f] = ions_.getCompileTimeResourcesViewList();
            d_c.setBuffer(&ionChargeDensity);
            d_m.setBuffer(&ionMassDensity);
            Vi.set_on(V);
            chargeFlux.set_on(c_f);
            ionTensor.set_on(m);
        }
        auto& populations = ions_.getRunTimeResourcesViewList();
        assert(populations.size() == nbrPops);
        for (std::size_t i = 0; i < nbrPops; ++i)
        {
            auto& buffers                       = popBuffers[i];
            auto const& [F, m, Np, Nc, poppack] = populations[i].getCompileTimeResourcesViewList();
            Np.setBuffer(&buffers.particleDensity);
            Nc.setBuffer(&buffers.chargeDensity);
            buffers.flux.set_on(F);
            buffers.momentumTensor.set_on(m);
            poppack.setBuffer(&buffers.pack);
        }
        return ions_;
    }


    TwoPopulationElectronsTest()
        : J{"J", layout, HybridQuantity::Vector::J}
        , Ve{"StandardHybridElectronFluxComputer_Ve", layout, HybridQuantity::Vector::V}
        , Vi{"bulkVel", layout, HybridQuantity::Vector::V}
        , chargeFlux{"chargeFlux", layout, HybridQuantity::Vector::V}
        , ionTensor{"momentumTensor", layout, HybridQuantity::Tensor::M}
        , ionChargeDensity{"chargeDensity", HybridQuantity::Scalar::rho,
                           layout.allocSize(HybridQuantity::Scalar::rho)}
        , ionMassDensity{"massDensity", HybridQuantity::Scalar::rho,
                         layout.allocSize(HybridQuantity::Scalar::rho)}
        , Pe{"Pe", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)}
        , popBuffers{{{pops[0].name, layout}, {pops[1].name, layout}}}
        , ions{makeIons()}
        , electrons{dict["electrons"], ions, J}
    {
        auto&& emm = std::get<0>(electrons.getCompileTimeResourcesViewList());
        auto&& fc  = std::get<0>(emm.getCompileTimeResourcesViewList());
        Ve.set_on(std::get<0>(fc.getCompileTimeResourcesViewList()));

        auto&& pc          = std::get<1>(emm.getCompileTimeResourcesViewList());
        auto const& [_, P] = pc.getCompileTimeResourcesViewList();
        P.setBuffer(&Pe);

        // uniform moments everywhere, ghost nodes included
        auto fill = [](auto& field, double const value) {
            std::fill(std::begin(field), std::end(field), value);
        };

        for (std::size_t i = 0; i < nbrPops; ++i)
        {
            auto const& pop = pops[i];
            auto& buffers   = popBuffers[i];
            fill(buffers.particleDensity, pop.density);
            for (std::size_t c = 0; c < 3; ++c)
                fill(buffers.flux[c], pop.density * pop.velocity[c]);
        }
        for (std::size_t c = 0; c < 3; ++c)
            fill(J[c], current[c]);
    }


    GridYee layout                     = NDlayout<dim, Base::interp>::create();
    PHARE::initializer::PHAREDict dict = createTwoPopulationDict(pops);

    UsableVecField<dim> J, Ve, Vi, chargeFlux;
    UsableTensorField<dim> ionTensor;
    GridND ionChargeDensity, ionMassDensity, Pe;
    std::array<PopulationBuffers, nbrPops> popBuffers;

    IonsT ions;
    Electrons<IonsT> electrons;
};


TYPED_TEST_SUITE(TwoPopulationElectronsTest, ElectronsTupleInfos, );


TYPED_TEST(TwoPopulationElectronsTest, ThatElectronVelocityIsChargeWeighted)
{
    auto& ions      = this->ions;
    auto& electrons = this->electrons;
    auto& layout    = this->layout;

    ions.computeChargeDensity();
    ions.computeBulkVelocity();
    ions.computeChargeFlux();
    electrons.update(layout);

    double Ne = 0;
    std::array<double, 3> chargeFlux{0., 0., 0.};
    for (auto const& pop : this->pops)
    {
        Ne += pop.charge * pop.density;
        for (std::size_t c = 0; c < 3; ++c)
            chargeFlux[c] += pop.charge * pop.density * pop.velocity[c];
    }

    for (std::size_t c = 0; c < 3; ++c)
    {
        auto const expected = (chargeFlux[c] - this->current[c]) / Ne;
        auto& Vecomp        = this->Ve[c];
        layout.evalOnBox(Vecomp, [&](auto const&... args) {
            EXPECT_DOUBLE_EQ(Vecomp(args...), expected) << "component " << c;
        });
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
