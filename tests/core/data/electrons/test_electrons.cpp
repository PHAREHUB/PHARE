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

const double Te = 0.12;



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

    dict["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    dict["electrons"]["pressure_closure"]["Te"]   = Te;

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
struct ElectronsTest : public ::testing::Test
{
    static constexpr auto dim         = typename TypeInfo::first_type{}();
    static constexpr auto interp      = typename TypeInfo::second_type{}();
    static constexpr auto densityName = std::string_view{"rho"};


    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;

    using GridND           = Grid<NdArrayVector<dim, floater_t<4>>, HybridQuantity::Scalar>;
    using FieldND          = Field<dim, HybridQuantity::Scalar, floater_t<4>>;
    using VecFieldND       = VecField<FieldND, HybridQuantity>;
    using SymTensorFieldND = SymTensorField<FieldND, HybridQuantity>;
    using ParticleArray_t  = ParticleArray<dim>;
    using IonPopulationND  = IonPopulation<ParticleArray_t, VecFieldND, SymTensorFieldND>;
    using IonsT            = Ions<IonPopulationND, GridYee>;
    using PartPackND       = ParticlesPack<ParticleArray_t>;
    using StandardHybridElectronFluxComputerT = StandardHybridElectronFluxComputer<IonsT>;


    GridYee layout = NDlayout<dim, interp>::create();

    Electromag<VecFieldND> electromag;

    UsableVecField<dim> J, F, Ve, Vi;
    UsableTensorField<dim> M, protons_M;

    GridND Nibuffer, NiProtons, Pe;

    ParticleArray_t domainParticles{layout.AMRBox()};
    ParticleArray_t patchGhostParticles = domainParticles;
    ParticleArray_t levelGhostParticles = domainParticles;
    PartPackND pack{"particles", &domainParticles, &patchGhostParticles, &levelGhostParticles};

    IonsT ions;
    Electrons<IonsT> electrons;

    template<typename... Args>
    auto static _ions(Args&... args)
    {
        auto const& [Fi, Nibuffer, NiProtons, Vi, M, protons_M, pack]
            = std::forward_as_tuple(args...);
        IonsT ions{createDict<dim>()["ions"]};
        {
            auto const& [V, m, d, md] = ions.getCompileTimeResourcesViewList();
            d.setBuffer(&Nibuffer);
            Vi.set_on(V);
            M.set_on(m);
        }
        auto& pops = ions.getRunTimeResourcesViewList();
        assert(pops.size() == 1);

        auto const& [F, m, d, poppack] = pops[0].getCompileTimeResourcesViewList();
        d.setBuffer(&NiProtons);
        Fi.set_on(F);
        protons_M.set_on(m);
        poppack.setBuffer(&pack);
        return ions;
    }


    ElectronsTest()
        : electromag{createDict<dim>()["electromag"]}
        , J{"J", layout, HybridQuantity::Vector::J}
        , F{"protons_flux", layout, HybridQuantity::Vector::V}
        , Ve{"StandardHybridElectronFluxComputer_Ve", layout, HybridQuantity::Vector::V}
        , Vi{"bulkVel", layout, HybridQuantity::Vector::V}
        , M{"momentumTensor", layout, HybridQuantity::Tensor::M}
        , protons_M{"protons_momentumTensor", layout, HybridQuantity::Tensor::M}
        , Nibuffer{std::string{densityName}, HybridQuantity::Scalar::rho,
                   layout.allocSize(HybridQuantity::Scalar::rho)}
        , NiProtons{"protons_rho", HybridQuantity::Scalar::rho,
                    layout.allocSize(HybridQuantity::Scalar::rho)}
        , Pe{"Pe", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)}
        , ions{_ions(F, Nibuffer, NiProtons, Vi, M, protons_M, pack)}
        , electrons{createDict<dim>()["electrons"], ions, J}
    {
        auto&& emm = std::get<0>(electrons.getCompileTimeResourcesViewList());
        auto&& fc  = std::get<0>(emm.getCompileTimeResourcesViewList());


        Ve.set_on(std::get<0>(fc.getCompileTimeResourcesViewList()));


        auto&& pc          = std::get<1>(emm.getCompileTimeResourcesViewList());
        auto const& [_, P] = pc.getCompileTimeResourcesViewList();
        P.setBuffer(&Pe);

        auto const& [Jx, Jy, Jz]    = J();
        auto const& [Vix, Viy, Viz] = Vi();


        if constexpr (dim == 1)
        {
            auto fill = [this](auto& field, auto const& filler) {
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
            auto fill = [this](auto& field, auto const& filler) {
                auto gsi_X = this->layout.ghostStartIndex(field, Direction::X);
                auto gei_X = this->layout.ghostEndIndex(field, Direction::X);
                auto gsi_Y = this->layout.ghostStartIndex(field, Direction::Y);
                auto gei_Y = this->layout.ghostEndIndex(field, Direction::Y);

                for (auto ix = gsi_X; ix <= gei_X; ++ix)
                {
                    for (auto iy = gsi_Y; iy <= gei_Y; ++iy)
                    {
                        auto point = this->layout.fieldNodeCoordinates(
                            field, Point<double, 2>{0., 0.}, ix, iy);
                        field(ix, iy) = filler(point[0], point[1]);
                    }
                }
            };

            fill(Vix, [](double x, double y) { return std::cosh(0.2 * x) * std::cosh(0.2 * y); });
            fill(Viy, [](double x, double y) { return std::cosh(0.3 * x) * std::cosh(0.3 * y); });
            fill(Viz, [](double x, double y) { return std::cosh(0.4 * x) * std::cosh(0.4 * y); });

            fill(Jx, [](double x, double y) { return std::sinh(0.2 * x) * std::sinh(0.2 * y); });
            fill(Jy, [](double x, double y) { return std::sinh(0.3 * x) * std::sinh(0.3 * y); });
            fill(Jz, [](double x, double y) { return std::sinh(0.4 * x) * std::sinh(0.4 * y); });

            fill(Nibuffer,
                 [](double x, double y) { return std::cosh(0.1 * x) * std::cosh(0.1 * y); });
        }
        else if constexpr (dim == 3)
        {
            auto fill = [this](auto& field, auto const& filler) {
                auto gsi_X = this->layout.ghostStartIndex(field, Direction::X);
                auto gei_X = this->layout.ghostEndIndex(field, Direction::X);
                auto gsi_Y = this->layout.ghostStartIndex(field, Direction::Y);
                auto gei_Y = this->layout.ghostEndIndex(field, Direction::Y);
                auto gsi_Z = this->layout.ghostStartIndex(field, Direction::Z);
                auto gei_Z = this->layout.ghostEndIndex(field, Direction::Z);

                for (auto ix = gsi_X; ix <= gei_X; ++ix)
                {
                    for (auto iy = gsi_Y; iy <= gei_Y; ++iy)
                    {
                        for (auto iz = gsi_Z; iz <= gei_Z; ++iz)
                        {
                            auto point = this->layout.fieldNodeCoordinates(
                                field, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                            field(ix, iy, iz) = filler(point[0], point[1], point[2]);
                        }
                    }
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

            fill(Nibuffer, [](double x, double y, double z) {
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
    auto& Ni = ions.density();

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

    electrons.update(layout);

    auto& Ne = electrons.density();

    auto check = [&layout](auto const& Vecomp, auto const& Vicomp, auto const& Jcomp,
                           auto const& Ne_, auto const& projector) -> void {
        if constexpr (dim == 1)
        {
            auto psi_X = layout.physicalStartIndex(Vicomp, Direction::X);
            auto pei_X = layout.physicalEndIndex(Vicomp, Direction::X);

            for (std::uint32_t i = psi_X; i < pei_X; ++i)
            {
                auto const JOnV = GridYee::project(Jcomp, {i}, projector());

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
            {
                for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
                {
                    auto const JOnV = GridYee::project(Jcomp, {i, j}, projector());

                    EXPECT_DOUBLE_EQ(Vecomp(i, j), Vicomp(i, j) - JOnV / Ne_(i, j));
                }
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
            {
                for (std::uint32_t j = psi_Y; j < pei_Y; ++j)
                {
                    for (std::uint32_t k = psi_Z; k < pei_Z; ++k)
                    {
                        auto const JOnV = GridYee::project(Jcomp, {i, j, k}, projector());

                        EXPECT_DOUBLE_EQ(Vecomp(i, j, k), Vicomp(i, j, k) - JOnV / Ne_(i, j, k));
                    }
                }
            }
        }
    };

    auto const& [Jx, Jy, Jz]    = this->J();
    auto const& [Vix, Viy, Viz] = this->Vi();
    auto const& [Vex, Vey, Vez] = this->Ve();
    check(Vex, Vix, Jx, Ne, &GridYee::JxToMoments);
    check(Vey, Viy, Jy, Ne, &GridYee::JyToMoments);
    check(Vez, Viz, Jz, Ne, &GridYee::JzToMoments);
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




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
