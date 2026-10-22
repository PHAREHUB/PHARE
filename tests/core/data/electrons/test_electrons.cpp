

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
#include <cmath>
#include <functional>
#include <complex>
#include <numbers>

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

    dict["electrons"]["pressure_closure"]["name"]  = std::string{"isothermal"};
    dict["electrons"]["pressure_closure"]["Te"]    = Te;
    dict["electrons"]["pressure_closure"]["Gamma"] = 1.;
    dict["electrons"]["pressure_closure"]["Pe"]    = static_cast<InitFunctionT>(ex); // todo?

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
struct ElectronsFixture

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

    UsableVecField<dim> B, J, F, Ve, Vi, chargeFlux;
    UsableTensorField<dim> ionTensor, protonTensor;

    GridND ionChargeDensity, ionMassDensity, protonParticleDensity, protonChargeDensity, Pe;

    GridND Te;


    ParticleArray_t domainParticles{layout.AMRBox()};
    ParticleArray_t patchGhostParticles = domainParticles;
    ParticleArray_t levelGhostParticles = domainParticles;
    PartPackND pack{"particles", &domainParticles, &patchGhostParticles, &levelGhostParticles};

    IonsT ions;
    StandardHybridElectronFluxComputerT fluxCompute;
    Electrons<StandardHybridElectronFluxComputerT> electrons;

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

    void initialize_variant_resources(auto& pressure_closure)
    {
        auto const visitors = varient_visitor_overloads{
            [&](VecFieldND& vf) {
                EXPECT_TRUE(vf.name() == "B");
                B.set_on(vf);
            },
            [&](FieldND& f) {
                if(f.name() == "Te") // could be tmp field
                    f.setBuffer(&Te);
            },
            [](auto&) {
                // if this happens you are missing an overload
                throw std::runtime_error("should not happen");
            },
        };
        for (auto& var : pressure_closure.getRunTimeResourcesViewList())
            std::visit(visitors, var);
    }


    ElectronsFixture(PHARE::initializer::PHAREDict const& dict = createDict<dim>())
        : electromag{dict["electromag"]}
        , B{"B", layout, HybridQuantity::Vector::B}
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
        , Te{"Te", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)}
        , ions{_ions(F, ionChargeDensity, ionMassDensity, protonParticleDensity,
                     protonChargeDensity, Vi, chargeFlux, ionTensor, protonTensor, pack)}
        , fluxCompute{ions, J}
        , electrons{dict["electrons"], fluxCompute, B}
    {
        /* TODO explain why... we have 2 flux computer : 1 is the flux computer and the same is a
         * copy in the pressure closure */
        auto&& emm = std::get<0>(electrons.getCompileTimeResourcesViewList());
        auto&& fc  = std::get<0>(emm.getCompileTimeResourcesViewList());
        auto&& pc  = std::get<1>(emm.getCompileTimeResourcesViewList());
        auto&& fc_ = std::get<0>(pc.getCompileTimeResourcesViewList());
        auto&& pe  = std::get<1>(pc.getCompileTimeResourcesViewList());

        Ve.set_on(std::get<0>(fc.getCompileTimeResourcesViewList()));
        Ve.set_on(std::get<0>(fc_.getCompileTimeResourcesViewList()));
        initialize_variant_resources(pc);
        pe.setBuffer(&Pe);
        EXPECT_TRUE(pc.isUsable());


        auto const& [Jx, Jy, Jz] = J();
        auto const& [Fx, Fy, Fz] = chargeFlux();

        if constexpr (dim == 1)
        {
            auto fill = [this](auto& field, auto const& filler) {
                for (auto const [amr_idx, lcl_idx] : layout.ghost_amr_lcl_idx(field))
                {
                    auto const point = this->layout.fieldNodeCoordinates(field, amr_idx);
                    field(lcl_idx)   = filler(point[0]);
                }
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
                for (auto const amr_idx : layout.AMRGhostBoxFor(field))
                {
                    auto point = this->layout.fieldNodeCoordinates(field, amr_idx);
                    field(layout.AMRToLocal(amr_idx)) = filler(point[0], point[1]);
                }
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
                for (auto const amr_idx : layout.AMRGhostBoxFor(field))
                {
                    auto point = this->layout.fieldNodeCoordinates(field, amr_idx);
                    field(layout.AMRToLocal(amr_idx)) = filler(point[0], point[1], point[2]);
                }
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


template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct ElectronsTest : public ElectronsFixture<TypeInfo>, public ::testing::Test
{
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
    auto const dt   = 0.0;

    electrons.update(layout, dt);

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
    auto const dt   = 0.0;

    electrons.update(layout, dt);


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
    auto const dt   = 0.0;

    electrons.update(layout, dt);

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


TEST(ElectronsFactoryTest, ThatThingsAreAsExpectedForCGL)
{
    auto dict = createDict<1>();

    dict["electrons"]["pressure_closure"]["name"] = std::string{"CGL"};

    ElectronsFixture<std::pair<DimConst<1>, InterpConst<1>>> fixture{dict};

    auto&& emm = std::get<0>(fixture.electrons.getCompileTimeResourcesViewList());
    auto&& pc  = std::get<1>(emm.getCompileTimeResourcesViewList());

    auto& B = pc.B();
}

TEST(ElectronsFactoryTest, ThatConstThingsAreAsExpectedForCGL)
{
    auto dict = createDict<1>();

    dict["electrons"]["pressure_closure"]["name"] = std::string{"CGL"};

    ElectronsFixture<std::pair<DimConst<1>, InterpConst<1>>> fixture{dict};

    auto&& emm     = std::get<0>(fixture.electrons.getCompileTimeResourcesViewList());
    auto const& pc = std::get<1>(emm.getCompileTimeResourcesViewList());

    auto& B = pc.B();
    EXPECT_TRUE(B.isUsable());
}


TEST(ElectronsFactoryTest, ThatConstThingsAreAsExpectedForPolytropic)
{
    auto dict = createDict<1>();

    dict["electrons"]["pressure_closure"]["name"] = std::string{"polytropic"};

    ElectronsFixture<std::pair<DimConst<1>, InterpConst<1>>> fixture{dict};

    auto&& emm     = std::get<0>(fixture.electrons.getCompileTimeResourcesViewList());
    auto const& pc = std::get<1>(emm.getCompileTimeResourcesViewList());

    auto& Te = pc.Te();
    EXPECT_TRUE(Te.isUsable());

    // the polytropic closure needs its step state ghosts filled: no one-call update
    EXPECT_THROW(fixture.electrons.update(fixture.layout, 0.1), std::runtime_error);
}


TEST(ElectronsFactoryTest, ThatThereIsNoB)
{
    auto const dict = createDict<1>();

    ElectronsFixture<std::pair<DimConst<1>, InterpConst<1>>> fixture{dict};

    auto&& emm     = std::get<0>(fixture.electrons.getCompileTimeResourcesViewList());
    auto const& pc = std::get<1>(emm.getCompileTimeResourcesViewList());

    EXPECT_ANY_THROW(auto& B = pc.B(););
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
        : B{"B", layout, HybridQuantity::Vector::B}
        , J{"J", layout, HybridQuantity::Vector::J}
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
        , electrons{dict["electrons"], StandardHybridElectronFluxComputer<IonsT>{ions, J}, B}
    {
        auto&& emm = std::get<0>(electrons.getCompileTimeResourcesViewList());
        auto&& fc  = std::get<0>(emm.getCompileTimeResourcesViewList());
        Ve.set_on(std::get<0>(fc.getCompileTimeResourcesViewList()));

        auto&& pc               = std::get<1>(emm.getCompileTimeResourcesViewList());
        auto const& [fc_, P]    = pc.getCompileTimeResourcesViewList();
        P.setBuffer(&Pe);
        // the pressure closure holds its own copy of the flux computer
        Ve.set_on(std::get<0>(fc_.getCompileTimeResourcesViewList()));

        // uniform moments everywhere, ghost nodes included
        auto fill = [](auto& field, double const value) {
            std::fill(std::begin(field), std::end(field), value);
        };

        for (std::size_t i = 0; i < nbrPops; ++i)
        {
            auto const& pop = pops[i];
            auto& buffers   = popBuffers[i];
            fill(buffers.particleDensity, pop.density);
            fill(buffers.chargeDensity, pop.charge * pop.density);
            for (std::size_t c = 0; c < 3; ++c)
                fill(buffers.flux[c], pop.density * pop.velocity[c]);
        }
        for (std::size_t c = 0; c < 3; ++c)
            fill(J[c], current[c]);
    }


    GridYee layout                     = NDlayout<dim, Base::interp>::create();
    PHARE::initializer::PHAREDict dict = createTwoPopulationDict(pops);

    UsableVecField<dim> B, J, Ve, Vi, chargeFlux;
    UsableTensorField<dim> ionTensor;
    GridND ionChargeDensity, ionMassDensity, Pe;
    std::array<PopulationBuffers, nbrPops> popBuffers;

    IonsT ions;
    Electrons<StandardHybridElectronFluxComputer<IonsT>> electrons;
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
    electrons.update(layout, /*dt=*/0.);

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




TYPED_TEST(TwoPopulationElectronsTest, ThatParticlesMustCarryTheirPopulationCharge)
{
    auto& sodium    = this->ions[1];
    auto& particles = this->popBuffers[1].domainParticles;

    typename TestFixture::ParticleArray_t::Particle_t particle;
    particle.weight = 1.;
    particle.charge = sodium.charge();
    particle.iCell  = this->layout.AMRBox().lower.toArray();
    particles.push_back(particle);
    EXPECT_TRUE(sodium.particlesHaveCharge());

    particles[0].charge = 2. * sodium.charge();
    EXPECT_FALSE(sodium.particlesHaveCharge());
    PHARE_DEBUG_DO(EXPECT_THROW(this->ions.computeChargeFlux(), std::runtime_error);)
}




// the tests below drive the closures alone, through a flux computer that only
// hands out a density and an electron velocity set by the test

namespace closure_test
{
template<std::size_t dim>
using GridLayout_t = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, 1}>::Hybrid::GridLayout_t;

template<std::size_t dim>
using Grid_t = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;

template<std::size_t dim>
struct FakeElectronFlux
{
    using GridLayout = GridLayout_t<dim>;
    using VecField   = typename UsableVecField<dim>::Super;
    using Field      = Field_t<dim>;

    NO_DISCARD bool isUsable() const { return true; }
    NO_DISCARD bool isSettable() const { return false; }

    NO_DISCARD Field const& density() const { return *N; }
    NO_DISCARD VecField const& velocity() const { return *V; }

    Field const* N;
    VecField const* V;
};

template<std::size_t dim>
struct PolytropicClosureFixture
{
    using Field1 = Field_t<dim>;
    using GridT  = Grid_t<dim>;

    double gamma = 5. / 3.;
    static constexpr double dx    = 0.1;
    static constexpr std::uint32_t nbrCells = dim == 1 ? 50 : 20;
    static constexpr double length          = dx * nbrCells;

    PHARE::initializer::PHAREDict makeDict() const
    {
        PHARE::initializer::PHAREDict dict;
        dict["pressure_closure"]["Gamma"] = gamma;
        dict["pressure_closure"]["Pe"]
            = static_cast<PHARE::initializer::InitFunction<dim>>(makeSharedPtr<dim>());
        return dict;
    }

    static auto makeLayout()
    {
        std::array<double, dim> meshSize, origin;
        std::array<std::uint32_t, dim> cells;
        meshSize.fill(dx);
        origin.fill(0.);
        cells.fill(nbrCells);
        return GridLayout_t<dim>{meshSize, cells, PHARE::core::Point<double, dim>{origin}};
    }

    explicit PolytropicClosureFixture(double const gamma_ = 5. / 3.)
        : gamma{gamma_}
    {
        Nf.setBuffer(&N);
        closure.pressure().setBuffer(&Pe);

        auto const visitors = varient_visitor_overloads{
            [&](Field1& f) {
                if (f.name() == "Te")
                    f.setBuffer(&Te);
                else if (f.name() == "Se_n")
                    f.setBuffer(&SeN);
                else
                    throw std::runtime_error("unexpected field resource " + f.name());
            },
            [](auto&) { throw std::runtime_error("unexpected resource"); },
        };
        for (auto& var : closure.getRunTimeResourcesViewList())
            std::visit(visitors, var);
    }

    GridLayout_t<dim> layout = makeLayout();

    // P, rho and V are all primal on Yee, so they share these indexes
    std::uint32_t const start = layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    std::uint32_t const end   = layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    // fn takes the node coordinates, ghost nodes included
    void fill(auto&& field, auto&& fn)
    {
        layout.evalOnGhostBox(field, [&](auto const&... ijk) {
            auto const index = std::array{ijk...};
            std::array<double, dim> x;
            for (std::size_t d = 0; d < dim; ++d)
                x[d] = (static_cast<double>(index[d]) - start) * dx;
            field(index) = std::apply(fn, x);
        });
    }

    // what the model's ghost filling schedule does on a periodic single patch domain:
    // a ghost node takes the domain node one period away, in every direction
    void periodicGhosts(GridT& grid)
    {
        auto const wrap = [&](std::uint32_t i) {
            if (i < start)
                return i + nbrCells;
            if (i > end)
                return i - nbrCells;
            return i;
        };
        layout.evalOnGhostBox(grid, [&](auto const&... ijk) {
            auto const index = std::array{ijk...};
            auto source      = index;
            for (auto& i : source)
                i = wrap(i);
            if (source != index)
                grid(index) = grid(source);
        });
    }

    // as HybridModel::fillElectronStepStateGhosts between beginStep and computePressure
    std::function<void(GridT&)> stepStateGhosts = [this](GridT& s) { periodicGhosts(s); };

    void step(double const dt, std::size_t const nbrStages = 1)
    {
        closure.beginStep(layout);
        stepStateGhosts(SeN);
        for (std::size_t stage = 0; stage < nbrStages; ++stage)
            closure.computePressure(layout, dt);
    }

    auto domainValues(GridT const& grid)
    {
        std::vector<double> values;
        layout.evalOnBox(grid, [&](auto const&... ijk) { values.push_back(grid(ijk...)); });
        return values;
    }

    GridT N{"N", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho), 1.};
    GridT Pe{"Pe", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P), 0.};
    GridT Te{"Te", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P), 0.};
    GridT SeN{"Se_n", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P), 0.};
    UsableVecField<dim> V{"Ve", layout, HybridQuantity::Vector::V};

    Field1 Nf{"N", HybridQuantity::Scalar::rho};
    PHARE::initializer::PHAREDict dict = makeDict();
    PolytropicElectronPressureClosure<FakeElectronFlux<dim>> closure{
        dict, FakeElectronFlux<dim>{&Nf, &V}};
};

using Fixture1D = PolytropicClosureFixture<1>;
using Fixture2D = PolytropicClosureFixture<2>;

} // namespace closure_test



// SolverPPC calls computePressure in predictor1, predictor2 and the corrector,
// each over the full dt: all must start from the beginStep state, or Pe
// would be advanced three times per time step
TEST(PolytropicElectronPressureClosure, stagesOfAStepDoNotAccumulate)
{
    closure_test::Fixture1D f;
    auto constexpr pi = std::numbers::pi;
    auto const Pe0    = [&](double x) { return 1. + 0.1 * std::sin(2 * pi * x / f.length); };

    f.fill(f.V(Component::X), [](double) { return 1.; });
    f.fill(f.Pe, Pe0);
    f.step(0.05);
    auto const oneStage = f.domainValues(f.Pe);

    f.fill(f.Pe, Pe0);
    f.step(0.05, /*nbrStages=*/3);
    auto const threeStages = f.domainValues(f.Pe);

    ASSERT_EQ(oneStage.size(), threeStages.size());
    for (std::size_t i = 0; i < oneStage.size(); ++i)
        EXPECT_DOUBLE_EQ(threeStages[i], oneStage[i]) << "at domain node " << i;

    // and the stage did move Pe, so the equality above is not trivially true
    f.fill(f.Pe, Pe0);
    EXPECT_NE(oneStage, f.domainValues(f.Pe));
}



// at Courant number 1 first order upwind is exact: s moves by exactly one node,
// downstream. An asymmetric profile makes a wrong upwind side visible.
TEST(PolytropicElectronPressureClosure, courantOneShiftsByExactlyOneNodeDownstream)
{
    auto constexpr pi = std::numbers::pi;

    for (double const v : {1., -1.})
    {
        SCOPED_TRACE("Ve = " + std::to_string(v));
        closure_test::Fixture1D f;
        auto const Pe0 = [&](double x) {
            return 1. + 0.1 * std::sin(2 * pi * x / f.length) + 0.05 * std::cos(4 * pi * x / f.length);
        };

        f.fill(f.V(Component::X), [&](double) { return v; });
        f.fill(f.Pe, Pe0);
        f.step(f.dx / std::abs(v));

        for (auto i = f.start; i <= f.end; ++i)
        {
            double const xi = (static_cast<double>(i) - f.start) * f.dx;
            EXPECT_NEAR(f.Pe(i), Pe0(xi - v / std::abs(v) * f.dx), 1e-13) << "at node " << i;
        }
    }
}



// for a single Fourier mode on a periodic grid the scheme has a closed form
// solution: s_i^n = 1 + a*Im(g^n exp(i k x_i)), with g the amplification factor
// of the upwind scheme. This checks transport speed, direction and damping.
TEST(PolytropicElectronPressureClosure, advectionMatchesTheExactDiscreteSolution)
{
    using namespace std::complex_literals;
    auto constexpr pi = std::numbers::pi;

    for (double const v : {1., -1.})
    {
        SCOPED_TRACE("Ve = " + std::to_string(v));
        closure_test::Fixture1D f;

        double const amplitude  = 0.1;
        double const wavenumber = 2 * pi * 5 / f.length; // 10 nodes per wavelength
        double const courant    = 0.5;
        double const dt         = courant * f.dx / std::abs(v);
        double const theta      = wavenumber * f.dx;
        std::size_t const nbrSteps = 200;

        // v > 0: s_i - c (s_i - s_{i-1});  v < 0: s_i + c (s_{i+1} - s_i)
        auto const g = v > 0 ? 1. - courant * (1. - std::exp(-1i * theta))
                             : 1. + courant * (std::exp(1i * theta) - 1.);

        f.fill(f.V(Component::X), [&](double) { return v; });
        f.fill(f.Pe, [&](double x) { return 1. + amplitude * std::sin(wavenumber * x); });

        for (std::size_t n = 0; n < nbrSteps; ++n)
            f.step(dt);

        auto const gn = std::pow(g, static_cast<double>(nbrSteps));
        for (auto i = f.start; i <= f.end; ++i)
        {
            double const xi       = (static_cast<double>(i) - f.start) * f.dx;
            double const expected = 1. + amplitude * std::imag(gn * std::exp(1i * wavenumber * xi));
            EXPECT_NEAR(f.Pe(i), expected, 1e-12) << "at node " << i;
        }
    }
}



// the unsplit update is a convex combination only if the Courant numbers summed
// over directions stay below 1: beyond, the closure must refuse to advance, and
// leave Pe as it was
TEST(PolytropicElectronPressureClosure, throwsBeyondTheAdvectionStabilityLimit)
{
    closure_test::Fixture1D f;
    f.fill(f.V(Component::X), [](double) { return 1.; });
    f.fill(f.Pe, [](double) { return 1.; });

    EXPECT_NO_THROW(f.step(1.0 * f.dx));

    auto const before = f.domainValues(f.Pe);
    EXPECT_THROW(f.step(1.5 * f.dx), std::runtime_error);
    EXPECT_EQ(before, f.domainValues(f.Pe));
}

TEST(PolytropicElectronPressureClosure, stabilityLimitSumsTheDirectionsIn2D)
{
    closure_test::Fixture2D f;
    f.fill(f.V(Component::X), [](double, double) { return 1.; });
    f.fill(f.V(Component::Y), [](double, double) { return 1.; });
    f.fill(f.Pe, [](double, double) { return 1.; });

    EXPECT_NO_THROW(f.step(0.5 * f.dx));  // 0.5 + 0.5
    EXPECT_THROW(f.step(0.75 * f.dx), std::runtime_error); // each 0.75 < 1, sum 1.5
}



// the 2D update is unsplit: for a single Fourier mode exp(i(kx x + ky y)) the
// amplification factor is 1 - sum_d of the 1D upwind terms. A flow with vx > 0 and
// vy < 0 checks both upwind sides and their combination.
TEST(PolytropicElectronPressureClosure, advectionMatchesTheExactDiscreteSolutionIn2D)
{
    using namespace std::complex_literals;
    auto constexpr pi = std::numbers::pi;

    closure_test::Fixture2D f;

    double const amplitude = 0.1;
    double const vx = 1., vy = -0.5;
    double const dt          = 0.4 * f.dx / (std::abs(vx) + std::abs(vy)); // sum of Courant 0.4
    double const kx          = 2 * pi * 2 / f.length;
    double const ky          = 2 * pi * 3 / f.length;
    double const cx          = std::abs(vx) * dt / f.dx;
    double const cy          = std::abs(vy) * dt / f.dx;
    std::size_t const nbrSteps = 100;

    auto const g = 1. - cx * (1. - std::exp(-1i * kx * f.dx)) + cy * (std::exp(1i * ky * f.dx) - 1.);

    f.fill(f.V(Component::X), [&](double, double) { return vx; });
    f.fill(f.V(Component::Y), [&](double, double) { return vy; });
    f.fill(f.Pe, [&](double x, double y) { return 1. + amplitude * std::sin(kx * x + ky * y); });

    for (std::size_t n = 0; n < nbrSteps; ++n)
        f.step(dt);

    auto const gn = std::pow(g, static_cast<double>(nbrSteps));
    f.layout.evalOnBox(f.Pe, [&](auto const&... ijk) {
        auto const index = std::array{ijk...};
        double const x   = (static_cast<double>(index[0]) - f.start) * f.dx;
        double const y   = (static_cast<double>(index[1]) - f.start) * f.dx;
        double const expected
            = 1. + amplitude * std::imag(gn * std::exp(1i * (kx * x + ky * y)));
        EXPECT_NEAR(f.Pe(index), expected, 1e-12) << "at node " << index[0] << "," << index[1];
    });
}



// with a uniform s (Pe = P0 (N/N0)^gamma) the transport is inert and the closure
// reduces to the algebraic polytropic law: Pe/N^gamma stays P0 whatever N does.
// Here V = a*x and N follows dN/dt = -a*N exactly.
TEST(PolytropicElectronPressureClosure, uniformEntropyGivesThePolytropicLaw)
{
    closure_test::Fixture1D f;

    double const a     = 0.1;
    double const dt    = 0.01;
    double const P0    = 1.;
    auto const density = [&](double t) { return std::exp(-a * t); };

    f.fill(f.V(Component::X), [&](double x) { return a * x; });
    f.fill(f.Pe, [&](double) { return P0; });

    // V = a*x is not periodic: the uniform state's ghosts are its nearest domain value
    auto const extrapolate = [&](auto& grid) {
        for (std::size_t i = 0; i < f.start; ++i)
            grid(i) = grid(f.start);
        for (std::size_t i = f.end + 1; i < grid.size(); ++i)
            grid(i) = grid(f.end);
    };

    double t = 0;
    for (std::size_t n = 0; n < 100; ++n)
    {
        // as in SolverPPC: moments at time n for beginStep, at n+1 for the corrector
        f.fill(f.N, [&](double) { return density(t); });
        f.closure.beginStep(f.layout);
        extrapolate(f.SeN);

        t += dt;
        f.fill(f.N, [&](double) { return density(t); });
        f.closure.computePressure(f.layout, dt);
    }

    for (auto i = f.start; i <= f.end; ++i)
        EXPECT_NEAR(f.Pe(i), P0 * std::pow(density(t), f.gamma), 1e-12) << "at node " << i;
}



// the closure must use the Gamma it is given, not a hard coded 5/3
TEST(PolytropicElectronPressureClosure, usesTheGivenGamma)
{
    closure_test::Fixture1D f{1.4};

    f.fill(f.V(Component::X), [](double) { return 0.; });
    f.fill(f.Pe, [](double) { return 2.; });
    f.closure.beginStep(f.layout);
    f.stepStateGhosts(f.SeN);

    f.fill(f.N, [](double) { return 0.5; }); // density halved within the step
    f.closure.computePressure(f.layout, 0.1);

    for (auto i = f.start; i <= f.end; ++i)
        EXPECT_NEAR(f.Pe(i), 2. * std::pow(0.5, 1.4), 1e-12) << "at node " << i;
}



// CGL is not implemented: it must refuse to run rather than produce a pressure
TEST(CGLElectronPressureClosure, refusesToRun)
{
    using Fixture = closure_test::Fixture1D;
    Fixture f;
    UsableVecField<1> B{"B", f.layout, HybridQuantity::Vector::B};

    CGLElectronPressureClosure<closure_test::FakeElectronFlux<1>> cgl{
        f.dict, closure_test::FakeElectronFlux<1>{&f.Nf, &f.V}, B};

    EXPECT_THROW(cgl.initialize(f.layout), std::runtime_error);
    EXPECT_THROW(cgl.computePressure(f.layout, 0.1), std::runtime_error);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
