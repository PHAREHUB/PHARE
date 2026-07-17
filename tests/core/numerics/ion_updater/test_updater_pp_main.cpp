// tests/core/numerics/ion_updater/test_updater_pp_main.cpp
//
//  requires
//  - cmake:    -DwithPhlop
//  - cxxflags: -DPHARE_LOG_LEVEL=1
//  - env:      PHARE_SCOPE_TIMING=1

#include "core/data/particles/particle_array_sorter.hpp"
#include "core/def/detail/mkn_avx.hpp"
#include <memory>
#define PHARE_SKIP_MPI_IN_CORE
#define PHARE_UNDEF_ASSERT
// #define PHARE_LOG_LEVEL 3
// #define PHARE_LOG_SCOPE_PRINT 1


#include "core/utilities/types.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/numerics/ion_updater/ion_updater_pc.hpp"
#include "core/numerics/ion_updater/ion_updater_per_particle.hpp"

#include "core/data/particles/particle_array_serializer.hpp"

#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "gtest/gtest.h"

namespace PHARE::core
{

template<auto opts, typename internals> // used by gtest
void PrintTo(ParticleArray<opts, internals> const& arr, std::ostream* os)
{
    // assert(arr.size());
    *os << arr;
}

bool constexpr static WITH_PATCH_GHOST  = false;
bool constexpr static USE_SERIALIZATION = 0;

auto static const bytes   = get_env_as("PHARE_GPU_BYTES", std::uint64_t{500000000});
auto static const cells   = get_env_as("PHARE_CELLS", std::uint32_t{4});
auto static const ppc     = get_env_as("PHARE_PPC", std::size_t{4});
auto static const seed    = get_env_as("PHARE_SEED", std::size_t{1012});
auto static const dt      = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const shufle  = get_env_as("PHARE_UNSORTED", std::size_t{0});
auto static const do_cmp  = get_env_as("PHARE_COMPARE", std::size_t{!USE_SERIALIZATION});
auto static const gen_bin = get_env_as("PHARE_GENERATE_PARTICLES_BIN", std::size_t{0});

std::string static const aos_particles_bin = "aos_particles.bin";

bool static const premain = []() {
    PHARE_WITH_MKN_GPU(                          //
        mkn::gpu::setLimitMallocHeapSize(bytes); //
    )
    PHARE_WITH_PHLOP( //
        PHARE_LOG_LINE_STR("cells: " << cells); PHARE_LOG_LINE_STR("ppc  : " << ppc);
        PHARE_LOG_LINE_STR("seed : " << seed);

        using namespace PHARE; //
        using namespace std::literals;
        if (auto e = core::get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
            phlop::threaded::ScopeTimerMan::INSTANCE()
                .file_name(".phare_times.0.txt")
                // .force_strings()
                // .headers("fn"s, "dim"s, "layout"s, "alloc"s, "storage"s, "time"s)
                .init(); //
    )
    return true;
}();

template<typename IonUpdater_t, typename Ions, typename EM, typename GridLayout_t>
auto construct_(Ions& ions, EM const& em, GridLayout_t const& layout)
{
    PHARE::initializer::PHAREDict dict;
    dict["simulation"]["algo"]["ion_updater"]["pusher"]["name"] = std::string{"modified_boris"};
    return IonUpdater_t{dict["simulation"]["algo"]["ion_updater"]};
}

template<typename Ions, typename EM, typename GridLayout_t>
auto get_updater_for(Ions& ions, EM const& em, GridLayout_t const& layout)
{
    using Particles = typename Ions::particle_array_type;
    if constexpr (any_in(Particles::layout_mode, LayoutMode::AoSMapped))
        return construct_<IonUpdater<Ions, EM, GridLayout_t>>(ions, em, layout);
    else if constexpr (Particles::layout_mode == LayoutMode::AoSPC)
        return construct_<IonUpdaterPC<Ions, EM, GridLayout_t>>(ions, em, layout);
    else
        return construct_<IonUpdaterPP<Ions, EM, GridLayout_t>>(ions, em, layout);
}


namespace detail::strings
{
    constexpr static std::string_view update = "update,";
    constexpr static std::string_view cma    = ",";
} // namespace detail::strings


template<typename Ions, typename EM, typename GridLayout_t>
void udpate(Ions& ions, EM const& em, GridLayout_t const& layout)
{
    // no timings
    auto updater       = get_updater_for(ions, em, layout);
    using IonUpdater_t = std::decay_t<decltype(updater)>;
    using Boxing_t = std::decay_t<decltype(*selection_boxing_impl<IonUpdater_t, GridLayout_t>())>;

    Boxing_t const boxing{layout, grow(layout.AMRBox(), GridLayout_t::nbrParticleGhosts())};
    updater.updatePopulations(ions, em, boxing, dt);
}

template<typename Ions, typename EM, typename GridLayout_t>
void update(Ions& ions, EM const& em, GridLayout_t const& layout)
{
    using Particles = typename Ions::particle_array_type;
    auto constexpr function_id
        = join_string_views_v<detail::strings::update, Particles::type_id, detail::strings::cma>;
    PHARE_LOG_LINE_STR(function_id);
    PHARE_LOG_SCOPE(1, function_id);


    udpate(ions, em, layout);
}


template<typename Particles_t, typename GridLayout_t>
auto make_ions(GridLayout_t const& layout)
{
    auto constexpr static alloc_mode = Particles_t::alloc_mode;
    auto constexpr static dim        = GridLayout_t::dimension;
    auto constexpr static interp     = GridLayout_t::interp_order;

    auto ions_p = std::make_shared<UsableIons_t<Particles_t, interp>>(layout, "protons");
    auto& ions  = *ions_p;

    if constexpr (USE_SERIALIZATION)
        if (!gen_bin)
            return ions_p; // NO REF!

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0);

    auto disperse = [](auto& particles) {
        delta_disperse(particles.domain_particles, seed);
        if constexpr (WITH_PATCH_GHOST)
            delta_disperse(particles.patch_ghost_particles, seed);
        if (shufle > 0)
        {
            shuffle(particles.domain_particles, seed);
            if constexpr (WITH_PATCH_GHOST)
                shuffle(particles.patch_ghost_particles, seed);
        }
    };

    auto add_particles = [&](auto& particles) {
        add_particles_in(particles.domain_particles, layout.AMRBox(), ppc);
        if constexpr (WITH_PATCH_GHOST)
            add_ghost_particles(particles.patch_ghost_particles, layout.AMRBox(), ppc,
                                GridLayout_t::nbrParticleGhosts());
    };

    add_particles(ions.populations[0].particles);
    disperse(ions.populations[0].particles);

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), layout.AMRBox().size() * ppc);

    if constexpr (USE_SERIALIZATION)
        if (gen_bin)
            serialize_particles(aos_particles_bin, ions.populations[0].particles.domain_particles);

    return ions_p;
}



template<typename Particles_t, typename GridLayout_t, typename Ions>
auto from_ions(GridLayout_t const& layout, Ions const& from)
{
    using FromParticles_t            = typename Ions::particle_array_type;
    auto constexpr static alloc_mode = Particles_t::alloc_mode;
    auto constexpr static dim        = GridLayout_t::dimension;
    auto constexpr static interp     = GridLayout_t::interp_order;
    auto ions_p = std::make_shared<UsableIons_t<Particles_t, interp>>(layout, "protons");
    auto& ions  = *ions_p;
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0);

    auto _add_particles_from = [&]<auto type>(auto& src, auto& dst) {
        ParticleArrayService::reserve_ppc_in<type>(dst, ppc);
        append_particles<type>(src, dst);
    };

    if constexpr (USE_SERIALIZATION)
        deserialize_particles<Particles_t, FromParticles_t>(
            aos_particles_bin, ions.populations[0].particles.domain_particles);
    else
        _add_particles_from.template operator()<ParticleType::Domain>(
            from.populations[0].particles.domain_particles,
            ions.populations[0].particles.domain_particles);

    if constexpr (WITH_PATCH_GHOST)
        _add_particles_from.template operator()<ParticleType::Ghost>(
            from.populations[0].particles.patch_ghost_particles,
            ions.populations[0].particles.patch_ghost_particles);

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), layout.AMRBox().size() * ppc);
    return ions_p;
}

template<bool quiet = false, typename Ions, typename GridLayout_t>
auto& evolve(Ions& ions, GridLayout_t const& layout)
{
    using Particles_t                = typename Ions::particle_array_type;
    auto constexpr static alloc_mode = Particles_t::alloc_mode;
    auto constexpr static dim        = GridLayout_t::dimension;

    assert(ions.populations[0].particles.domain_particles.size() > 0);
    UsableElectromag<GridLayout_t, alloc_mode> em{layout};
    if constexpr (quiet)
        udpate(*ions, *em, layout);
    else
        update(*ions, *em, layout);
    assert(ions.populations[0].particles.domain_particles.size() > 0);
    return ions;
}


template<typename GridLayout_t, typename P0, typename P1>
void check_particles(GridLayout_t const& layout, P0& ref, P1& cmp, std::size_t ghosts = 0)
{
    using CPU_ref = ParticleArray<ParticleArrayOptions{P0::dimension, LayoutMode::AoS,
                                                       StorageMode::VECTOR, AllocatorMode::CPU}>;

    auto const box = grow(layout.AMRBox(), ghosts);

    sort_particles(cmp, box);
    sort_particles(ref, box);

    EXPECT_EQ(ref.size(), cmp.size());

    auto const report = compare_particles(ref, cmp);
    if (!report)
    {
        PHARE_LOG_LINE_STR("Comparing Particle Arrays: " << P0::id() << " vs " << P1::id());
        PHARE_LOG_LINE_STR("results: " << report.why());
    }
    auto const tmp = convert_particles_and_sort<CPU_ref>(cmp, layout);
    EXPECT_EQ(ref, tmp);
    EXPECT_TRUE(report);
}

template<typename GridLayout_t, typename R, typename C>
void compare(GridLayout_t const& layout, R& ref, C& cmp)
{
    if (!do_cmp)
        return;

    using ParticleArray_t = typename C::ParticleArray_t;

    PHARE_LOG_LINE_STR(ref.populations[0].particles.domain_particles.begin().copy());
    PHARE_LOG_LINE_STR(cmp.populations[0].particles.domain_particles.begin().copy());

    check_particles(layout, ref.populations[0].particles.domain_particles,
                    cmp.populations[0].particles.domain_particles);

    if constexpr (WITH_PATCH_GHOST)
        check_particles(layout, ref.populations[0].particles.patch_ghost_particles,
                        cmp.populations[0].particles.patch_ghost_particles);
    using enum LayoutMode;
    using enum AllocatorMode;

    double diff = 1e-13; // 0.0000000000000011102230
    if constexpr (ParticleArray_t::alloc_mode == GPU_UNIFIED)
        diff *= 1e2; // atomics no order guaranteed#
    // else if constexpr (any_in(ParticleArray_t::layout_mode, AoSTS, SoATS))
    //     diff *= 1e1; // p2m op order diff
    // else if constexpr (any_in(ParticleArray_t::layout_mode, SoAVX))
    //     diff *= 1e1; // not sure really


    {
        auto const freport
            = compare_tensor_fields(*ref.populations[0].F, *cmp.populations[0].F, diff);
        // if (!freport)
        PHARE_LOG_LINE_STR("results: " << freport.why());
        EXPECT_TRUE(freport);
    }
    auto const freport = compare_fields(*ref.populations[0].rho, *cmp.populations[0].rho, diff);
    // if (!freport)
    PHARE_LOG_LINE_STR("results: " << freport.why());
    EXPECT_TRUE(freport);
}


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode = PHARE::AllocatorMode::CPU>
struct TestParam
{
    static_assert(std::is_same_v<decltype(_layout_mode), LayoutMode>);
    static_assert(std::is_same_v<decltype(_alloc_mode), PHARE::AllocatorMode>);
    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = _layout_mode;
    auto constexpr static alloc_mode  = _alloc_mode;
};

template<typename Ions, std::size_t dim, std::size_t interp>
struct DefaultIons
{
    using GridLayout_t = TestGridLayout<typename PHARE_Types<SimOpts{dim, interp}>::GridLayout_t>;

    GridLayout_t const layout{cells};
    std::shared_ptr<Ions> init; // = make_ions<typename Ions::particle_array_type>(layout);
    std::shared_ptr<Ions> ions; // = make_ions<typename Ions::particle_array_type>(layout);

    DefaultIons()
    {
        assert(premain);
        init = make_ions<typename Ions::particle_array_type>(layout);

        PHARE_LOG_LINE_SS("");
        if constexpr (!USE_SERIALIZATION)
        {
            PHARE_LOG_LINE_SS("");
            ions = make_ions<typename Ions::particle_array_type>(layout);
            evolve<true>(*ions, *layout);
        }
    }
};

template<std::size_t dim, std::size_t interp>
using DefaultIons_t = DefaultIons<UsableIons_t<AoSMappedParticleArray<dim>, interp>, dim, interp>;


// static inline auto refs = std::tuple<DefaultIons_t<1, 1>, int, int>{}; // only 1d
// static inline auto refs = std::tuple<int, DefaultIons_t<2, 1>, int>{};// only 2d
static inline auto refs = std::tuple<int, int, DefaultIons_t<3, 1>>{}; // only 3d

// static inline auto refs
//     = std::tuple<DefaultIons_t<1, 1>, DefaultIons_t<2, 1>, DefaultIons_t<3, 1>>{}; // all d


template<typename Param>
struct IonUpdaterPPTest : public ::testing::Test
{
    auto constexpr static dim         = Param::dim;
    auto constexpr static interp      = 1;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;


    using GridLayout_t = TestGridLayout<typename PHARE_Types<SimOpts{dim, interp}>::GridLayout_t>;
    using RefParticleArray_t = AoSMappedParticleArray<dim>;
    using CmpParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    using RefIons_t = UsableIons_t<RefParticleArray_t, interp>;
    using CmpIons_t = UsableIons_t<CmpParticleArray_t, interp>;
    using DefIons   = DefaultIons<RefIons_t, dim, interp>;

    DefIons& ref = std::get<dim - 1>(refs);
    GridLayout_t const layout{cells};


    IonUpdaterPPTest() {}

    auto make_ions() const { return from_ions<CmpParticleArray_t>(layout, *ref.init); }
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !
//      TestParam<1, LayoutMode::AoS>

// PHARE_WITH_THRUST(
//     ,TestParam<1, LayoutMode::SoA>
//     ,TestParam<1, LayoutMode::SoAVX>
// )

//     ,TestParam<2, LayoutMode::AoS>

// PHARE_WITH_THRUST(
//     ,TestParam<2, LayoutMode::SoA>
//     ,TestParam<2, LayoutMode::SoAVX>
// )

     TestParam<3, LayoutMode::AoS>
    ,TestParam<3, LayoutMode::AoSMapped>

PHARE_WITH_MKN_AVX(
    // ,TestParam<3, LayoutMode::SoAVX>
)


PHARE_WITH_THRUST(
    // ,TestParam<3, LayoutMode::SoA>
    // ,TestParam<3, LayoutMode::SoAVX>
)


PHARE_WITH_MKN_GPU(
    // ,TestParam<1, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>   // 3
    // ,TestParam<1, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED> // 4
    // ,TestParam<1, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED,  /*impl=*/1> // 5
    // ,TestParam<1, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED,  /*impl=*/2>
    // ,TestParam<1, LayoutMode::SoA, AllocatorMode::GPU_UNIFIED>
)
    // ,TestParam<2, LayoutMode::AoS>
    // ,TestParam<2, LayoutMode::AoSPC>
PHARE_WITH_MKN_GPU(
    // ,TestParam<2, LayoutMode::SoA>
    // ,TestParam<2, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>
    // ,TestParam<2, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED>
    // ,TestParam<2, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, /*impl=*/1> // 13
    // ,TestParam<2, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, /*impl=*/2> // 14
    // ,TestParam<2, LayoutMode::SoA, AllocatorMode::GPU_UNIFIED>
)
    // /*,*/TestParam<3, LayoutMode::AoS>
    // ,TestParam<3, LayoutMode::AoSPC>
PHARE_WITH_THRUST(
    // ,TestParam<3, LayoutMode::SoA>
)

PHARE_WITH_MKN_GPU(
    ,TestParam<3, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>

//     // ,TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED>
//     // ,TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, /*impl=*/1> //
//     // ,TestParam<3, LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, /*impl=*/2> // 22
// PHARE_WITH_THRUST(
//     ,TestParam<3, LayoutMode::SoA, AllocatorMode::GPU_UNIFIED>
// )
)

>;
// clang-format on

TYPED_TEST_SUITE(IonUpdaterPPTest, Permutations_t);


TYPED_TEST(IonUpdaterPPTest, updater)
{
    if constexpr (PHARE::core::USE_SERIALIZATION)
        if (gen_bin)
            return;

    auto cmp_ions = this->make_ions();

    // if (do_cmp)
    // {
    //     PHARE_LOG_LINE_SS("Init check!")
    //     compare(*this->layout, *this->ref.init, *cmp_ions);
    // }

    compare(*this->layout,   //
            *this->ref.ions, //
            evolve(*cmp_ions, *this->layout));
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    auto r = RUN_ALL_TESTS();
    PHARE_WITH_PHLOP(phlop::threaded::ScopeTimerMan::reset());
    PHARE_WITH_MKN_AVX({
        for (auto const& [k, v] : mkn::avx::Counter::I().cnts)
        {
            PHARE_LOG_LINE_SS(k << " " << v);
        }
    })
    return r;
}
