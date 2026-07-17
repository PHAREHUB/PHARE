

// #include "core/hybrid/hybrid_quantities.hpp"
#include "phare_core.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "tests/initializer/init_functions.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

#include <cstdint>


using namespace PHARE;
using namespace PHARE::core;


std::size_t constexpr static interp    = 1;
std::size_t constexpr static cells     = 14;
std::size_t constexpr static tile_size = 4;


template<std::size_t _dim, auto lm, auto am>
struct TestParam
{
    static_assert(all_are<LayoutMode>(lm));
    static_assert(all_are<AllocatorMode>(am));

    auto constexpr static dim      = _dim;
    auto constexpr static sim_opts
        = SimOpts{.dimension = dim, .interp_order = interp, .layout_mode = lm, .alloc_mode = am};

    using PhareTypes   = core::PHARE_Types<sim_opts>;
    using Box_t        = PHARE::core::Box<int, dim>;
    using GridLayout_t = TestGridLayout<typename PhareTypes::GridLayout_t>;
};

struct InitFunctor
{
    using Param  = std::vector<double> const&;
    using Return = std::shared_ptr<PHARE::core::Span<double>>;

    Return operator()(Param x) { return std::make_shared<core::VectorSpan<double>>(x); }
    Return operator()(Param x, Param) { return std::make_shared<core::VectorSpan<double>>(x); }
    Return operator()(Param x, Param, Param)
    {
        return std::make_shared<core::VectorSpan<double>>(x);
    }

    std::size_t idx = 0;
};



template<std::size_t dim>
auto static em_dict()
{
    using InitFunctionT = PHARE::initializer::InitFunction<dim>;
    static InitFunctor bx, by, bz;

    initializer::PHAREDict dict;
    dict["name"]                                   = std::string{"EM"};
    dict["electric"]["name"]                       = std::string{"E"};
    dict["magnetic"]["name"]                       = std::string{"B"};
    dict["magnetic"]["initializer"]["x_component"] = static_cast<InitFunctionT>(bx);
    dict["magnetic"]["initializer"]["y_component"] = static_cast<InitFunctionT>(by);
    dict["magnetic"]["initializer"]["z_component"] = static_cast<InitFunctionT>(bz);

    return dict;
}


template<typename TestParam>
struct Patch
{
    auto constexpr static dim = TestParam::dim;

    using PhareTypes   = TestParam::PhareTypes;
    using GridLayout_t = TestParam::GridLayout_t;
    using Box_t        = TestParam::Box_t;

    using Field_t      = PhareTypes::Field_t;
    using Grid_t       = PhareTypes::Grid_t;
    using VecField_t   = PhareTypes::VecField_t;
    using Electromag_t = PhareTypes::Electromag_t;

    Patch(Box_t const& box)
        : layout{box}
    {
        for (std::size_t i = 0; i < 3; ++i)
        {
            em.E[i].setBuffer(&Es[i]);
            em.B[i].setBuffer(&Bs[i]);
        }
    }

    Grid_t make_grid(std::string const& name, auto qty) const { return {name, layout, qty}; }

    GridLayout_t const layout;

    std::array<Grid_t, 3> Bs{make_grid("EM_B_x", HybridQuantity::Scalar::Bx),
                             make_grid("EM_B_y", HybridQuantity::Scalar::By),
                             make_grid("EM_B_z", HybridQuantity::Scalar::Bz)};
    std::array<Grid_t, 3> Es{make_grid("EM_E_x", HybridQuantity::Scalar::Ex),
                             make_grid("EM_E_y", HybridQuantity::Scalar::Ey),
                             make_grid("EM_E_z", HybridQuantity::Scalar::Ez)};

    Electromag_t em{em_dict<dim>()};
};

template<std::size_t dim, typename TestParam>
struct AGridFieldTest;

template<typename TestParam>
struct AGridFieldTest<1, TestParam>
{
    using Box_t = TestParam::Box_t;

    AGridFieldTest()
    {
        patches.reserve(3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
        {
            auto const cellx = i * cells;
            patches.emplace_back(Box_t{Point{cellx}, Point{cellx + off}});
        }
    }

    std::vector<Patch<TestParam>> patches;
};


template<typename TestParam>
struct AGridFieldTest<2, TestParam>
{
    using Box_t = TestParam::Box_t;

    AGridFieldTest()
    {
        patches.reserve(3 * 3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
            {
                auto const cellx = i * cells;
                auto const celly = j * cells;
                patches.emplace_back(Box_t{Point{cellx, celly}, Point{cellx + off, celly + off}});
            }
    }

    std::vector<Patch<TestParam>> patches;
};


template<typename TestParam>
struct AGridFieldTest<3, TestParam>
{
    using Box_t = TestParam::Box_t;

    AGridFieldTest()
    {
        patches.reserve(3 * 3);
        auto const off = cells - 1;
        for (std::uint8_t i = 0; i < 3; ++i)
            for (std::uint8_t j = 0; j < 3; ++j)
                for (std::uint8_t k = 0; k < 3; ++k)
                {
                    auto const cellx = i * cells;
                    auto const celly = j * cells;
                    auto const cellz = k * cells;
                    patches.emplace_back(Box_t{Point{cellx, celly, cellz},
                                               Point{cellx + off, celly + off, cellz + off}});
                }
    }

    std::vector<Patch<TestParam>> patches;
};


template<typename TestParam_>
struct GridFieldTest : public ::testing::Test, public AGridFieldTest<TestParam_::dim, TestParam_>
{
    using TestParam    = TestParam_;
    using Super        = AGridFieldTest<TestParam::dim, TestParam>;
    using GridLayout_t = TestParam::GridLayout_t;

    using Super::patches;

    GridFieldTest() {}
};


// clang-format off
using ParticlesDatas = testing::Types< //

    TestParam<1, LayoutMode::AoS, AllocatorMode::CPU>

PHARE_WITH_MKN_GPU(
    ,TestParam<1, LayoutMode::AoSTS, AllocatorMode::CPU>
    ,TestParam<2, LayoutMode::AoSTS, AllocatorMode::CPU>
)

PHARE_WITH_GPU(
   ,TestParam<2, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>
)//PHARE_WITH_GPU

>;
// clang-format on

TYPED_TEST_SUITE(GridFieldTest, ParticlesDatas);


namespace PHARE::core
{

TYPED_TEST(GridFieldTest, test_compiles)
{
    for (auto& patch : this->patches)
        patch.em.initialize(patch.layout);
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
