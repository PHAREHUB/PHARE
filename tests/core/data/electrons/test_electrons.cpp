#include "initializer/data_provider.hpp"

#include "core/data/electrons/electrons.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/electromag/electromag.hpp"
#include "src/core/utilities/types.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "tests/initializer/init_functions.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"


#include <string_view>

using namespace PHARE::core;

double const Te = 0.12;



template<int dim>
PHARE::initializer::PHAREDict createDict()
{
    using InitFunctionT = PHARE::initializer::InitFunction<dim>;

    auto density = makeSharedPtr<dim>();
    auto vx      = makeSharedPtr<dim>();
    auto vy      = makeSharedPtr<dim>();
    auto vz      = makeSharedPtr<dim>();
    auto vthx    = makeSharedPtr<dim>();
    auto vthy    = makeSharedPtr<dim>();
    auto vthz    = makeSharedPtr<dim>();
    auto bx      = makeSharedPtr<dim>();
    auto by      = makeSharedPtr<dim>();
    auto bz      = makeSharedPtr<dim>();
    auto ex      = makeSharedPtr<dim>();
    auto ey      = makeSharedPtr<dim>();
    auto ez      = makeSharedPtr<dim>();


    PHARE::initializer::PHAREDict dict;
    dict["ions"]["nbrPopulations"]                       = std::size_t{1};
    dict["ions"]["pop0"]["name"]                         = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]                         = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["name"] = std::string{"maxwellian"};

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

    dict["ions"]["pop0"]["particle_initializer"]["nbr_part_per_cell"] = int{100};
    dict["ions"]["pop0"]["particle_initializer"]["charge"]            = -1.;
    dict["ions"]["pop0"]["particle_initializer"]["basis"]             = std::string{"cartesian"};

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

    using nDL = GridLayout<GridLayoutImplYee<dim, interp>>;

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


    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;

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

    UsableVecField<dim> B, J, F, Ve, Vi;
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
                     protonChargeDensity, Vi, ionTensor, protonTensor, pack]
            = std::forward_as_tuple(args...);
        IonsT ions{createDict<dim>()["ions"]};
        {
            auto const& [V, m, d_c, d_m] = ions.getCompileTimeResourcesViewList();
            d_c.setBuffer(&ionChargeDensity);
            d_m.setBuffer(&ionMassDensity);
            Vi.set_on(V);
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
                EXPECT_TRUE(f.name() == "Te");
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
                     protonChargeDensity, Vi, ionTensor, protonTensor, pack)}
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


        auto const& [Jx, Jy, Jz]    = J();
        auto const& [Vix, Viy, Viz] = Vi();

        if constexpr (dim == 1)
        {
            auto fill = [this](auto& field, auto const& filler) {
                for (auto const [amr_idx, lcl_idx] : layout.ghost_amr_lcl_idx(field))
                {
                    auto const point = this->layout.fieldNodeCoordinates(field, amr_idx);
                    field(lcl_idx)   = filler(point[0]);
                }
            };

            fill(Vix, [](double x) { return std::cosh(0.2 * x); });
            fill(Viy, [](double x) { return std::cosh(0.3 * x); });
            fill(Viz, [](double x) { return std::cosh(0.4 * x); });

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

            fill(Vix, [](double x, double y) { return std::cosh(0.2 * x) * std::cosh(0.2 * y); });
            fill(Viy, [](double x, double y) { return std::cosh(0.3 * x) * std::cosh(0.3 * y); });
            fill(Viz, [](double x, double y) { return std::cosh(0.4 * x) * std::cosh(0.4 * y); });

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

            fill(Vix, [](double x, double y, double z) {
                return std::cosh(0.2 * x) * std::cosh(0.2 * y) * std::cosh(0.2 * z);
            });
            fill(Viy, [](double x, double y, double z) {
                return std::cosh(0.3 * x) * std::cosh(0.3 * y) * std::cosh(0.3 * z);
            });
            fill(Viz, [](double x, double y, double z) {
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


TYPED_TEST(ElectronsTest, ThatElectronsVelocityEqualIonVelocityMinusJ)
{
    static constexpr auto dim    = typename TypeParam::first_type{}();
    static constexpr auto interp = typename TypeParam::second_type{}();
    using GridYee                = GridLayout<GridLayoutImplYee<dim, interp>>;

    auto& electrons = this->electrons;
    auto& layout    = this->layout;
    auto const dt   = 0.0;

    electrons.update(layout, dt);

    auto& Ne = electrons.density();

    auto check = [&layout]<auto projector>(auto const& Vecomp, auto const& Vicomp,
                                           auto const& Jcomp, auto const& Ne_) -> void {
        if constexpr (dim == 1)
        {
            auto psi_X = layout.physicalStartIndex(Vicomp, Direction::X);
            auto pei_X = layout.physicalEndIndex(Vicomp, Direction::X);
            for (std::uint32_t i = psi_X; i < pei_X; ++i)
            {
                auto const JOnV = GridYee::template project<projector>(Jcomp, {i});
                EXPECT_DOUBLE_EQ(Vecomp(i), Vicomp(i) - JOnV / Ne_(i));
            }
        }
        else if constexpr (dim == 2)
        {
            auto psi_X = layout.physicalStartIndex(Vicomp, Direction::X);
            auto pei_X = layout.physicalEndIndex(Vicomp, Direction::X);
            auto psi_Y = layout.physicalStartIndex(Vicomp, Direction::Y);
            auto pei_Y = layout.physicalEndIndex(Vicomp, Direction::Y);
            for (std::uint32_t i = psi_X; i < pei_X; ++i)
                for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
                {
                    auto const JOnV = GridYee::template project<projector>(Jcomp, {i, j});
                    EXPECT_DOUBLE_EQ(Vecomp(i, j), Vicomp(i, j) - JOnV / Ne_(i, j));
                }
        }
        else if constexpr (dim == 3)
        {
            auto psi_X = layout.physicalStartIndex(Vicomp, Direction::X);
            auto pei_X = layout.physicalEndIndex(Vicomp, Direction::X);
            auto psi_Y = layout.physicalStartIndex(Vicomp, Direction::Y);
            auto pei_Y = layout.physicalEndIndex(Vicomp, Direction::Y);
            auto psi_Z = layout.physicalStartIndex(Vicomp, Direction::Z);
            auto pei_Z = layout.physicalEndIndex(Vicomp, Direction::Z);
            for (std::uint32_t i = psi_X; i < pei_X; ++i)
                for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
                    for (std::uint32_t k = psi_Z; k < pei_Z; ++k)
                    {
                        auto const JOnV = GridYee::template project<projector>(Jcomp, {i, j, k});
                        EXPECT_DOUBLE_EQ(Vecomp(i, j, k), Vicomp(i, j, k) - JOnV / Ne_(i, j, k));
                    }
        }
    };

    auto const& [Jx, Jy, Jz]    = this->J();
    auto const& [Vix, Viy, Viz] = this->Vi();
    auto const& [Vex, Vey, Vez] = this->Ve();

    check.template operator()<GridYee::JxToMoments>(Vex, Vix, Jx, Ne);
    check.template operator()<GridYee::JyToMoments>(Vey, Viy, Jy, Ne);
    check.template operator()<GridYee::JzToMoments>(Vez, Viz, Jz, Ne);
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
}


TEST(ElectronsFactoryTest, ThatThereIsNoB)
{
    auto const dict = createDict<1>();

    ElectronsFixture<std::pair<DimConst<1>, InterpConst<1>>> fixture{dict};

    auto&& emm     = std::get<0>(fixture.electrons.getCompileTimeResourcesViewList());
    auto const& pc = std::get<1>(emm.getCompileTimeResourcesViewList());

    EXPECT_ANY_THROW(auto& B = pc.B(););
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
