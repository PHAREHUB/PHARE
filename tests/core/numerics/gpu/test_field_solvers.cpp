
#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"
#include "initializer/data_provider.hpp"

#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/data/electrons/electrons.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "gtest/gtest.h"

namespace PHARE::core
{
std::uint32_t static constexpr cells = 30;
double static constexpr dt           = .001;
double static constexpr eta          = .1;
double static constexpr nu           = .01;
double static constexpr Te           = .5;
double static constexpr atol         = 1e-14;

double max_abs(auto const& vecfield)
{
    double max = 0;
    for (std::size_t c = 0; c < 3; ++c)
        for (std::size_t i = 0; i < vecfield[c].size(); ++i)
            max = std::max(max, std::abs(vecfield[c].data()[i]));
    return max;
}

void fill(auto& field, double const offset)
{
    for (std::size_t i = 0; i < field.size(); ++i)
        field.data()[i] = 1 + offset + .01 * static_cast<double>(i % 7);
}

template<auto opts>
struct FieldSolvers
{
    using Hybrid                     = PHARE_Types<opts>::Hybrid;
    using GridLayout_t               = Hybrid::GridLayout_t;
    using IonTypes                   = UsableIonsDefaultTypes<opts>;
    using Ions_t                     = UsableIons<IonTypes>;
    using IonsBase_t                 = Ions<typename IonTypes::IonPopulation_t, GridLayout_t>;
    auto constexpr static field_opts = TensorFieldOptions<Hybrid>{};
    using Grid_t                     = decltype(field_opts)::Grid_t;

    static auto electrons_dict()
    {
        initializer::PHAREDict dict;
        dict["pressure_closure"]["name"] = std::string{"isothermal"};
        dict["pressure_closure"]["Te"]   = Te;
        return dict;
    }

    FieldSolvers(HyperMode const hyper_mode)
        : hyper_mode{hyper_mode}
    {
        fill(ions.chargeDensity(), 0);
        for (std::size_t c = 0; c < 3; ++c)
            fill(ions.chargeFlux()[c], .1 * static_cast<double>(c));

        auto&& emm = std::get<0>(electrons.getCompileTimeResourcesViewList());
        auto&& fc  = std::get<0>(emm.getCompileTimeResourcesViewList());
        Ve.set_on(std::get<0>(fc.getCompileTimeResourcesViewList()));
        auto const& [_, P]
            = std::get<1>(emm.getCompileTimeResourcesViewList()).getCompileTimeResourcesViewList();
        P.setBuffer(&Pe);
    }

    void run()
    {
        Faraday<GridLayout_t>{layout}(*em.B, *em.E, *emNew.B, dt);
        Ampere<GridLayout_t>{layout}(*emNew.B, *J);
        electrons.update(layout);
        Ohm<GridLayout_t>{{eta, nu, hyper_mode}, layout}(ions.chargeDensity(), *Ve, *Pe, *emNew.B,
                                                         *J, *emNew.E);
    }

    HyperMode const hyper_mode;
    GridLayout_t layout = TestGridLayout<GridLayout_t>{cells};
    UsableElectromag<field_opts> em{layout};
    UsableElectromag<field_opts> emNew{layout, false};
    UsableVecField<field_opts> J{"J", layout, HybridQuantity::Vector::J};
    UsableVecField<field_opts> Ve{"StandardHybridElectronFluxComputer_Ve", layout,
                                  HybridQuantity::Vector::V};
    Grid_t Pe{"Pe", layout, HybridQuantity::Scalar::P};
    Ions_t ions{layout, "protons"};
    Electrons<IonsBase_t> electrons{electrons_dict(), ions, J};
};


template<std::size_t dim, std::size_t interp>
void test(HyperMode const hyper_mode)
{
    auto constexpr static cpu_opts
        = SimOpts{.dimension = dim, .interp_order = interp, .alloc_mode = AllocatorMode::CPU};
    auto constexpr static gpu_opts = SimOpts{
        .dimension = dim, .interp_order = interp, .alloc_mode = AllocatorMode::GPU_UNIFIED};

    FieldSolvers<cpu_opts> cpu{hyper_mode};
    FieldSolvers<gpu_opts> gpu{hyper_mode};
    cpu.run();
    gpu.run();

    EXPECT_TRUE(cpu.emNew.B.isclose(gpu.emNew.B, atol));
    EXPECT_TRUE(cpu.J.isclose(gpu.J, atol));
    EXPECT_TRUE(cpu.Ve.isclose(gpu.Ve, atol));
    EXPECT_TRUE(cpu.Pe.isclose(gpu.Pe, atol));
    EXPECT_TRUE(cpu.emNew.E.isclose(gpu.emNew.E, atol * max_abs(cpu.emNew.E)));
}

TEST(FieldSolversGPU, matchesCPU_1d)
{
    test<1, 1>(HyperMode::constant);
    test<1, 3>(HyperMode::spatial);
}

TEST(FieldSolversGPU, matchesCPU_2d)
{
    test<2, 1>(HyperMode::constant);
    test<2, 3>(HyperMode::spatial);
}

TEST(FieldSolversGPU, matchesCPU_3d)
{
    test<3, 1>(HyperMode::constant);
    test<3, 3>(HyperMode::spatial);
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    PHARE_WITH_KOKKOS(Kokkos::initialize(argc, argv);)
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
