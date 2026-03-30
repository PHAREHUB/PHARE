

#include "core/data/particles/particle_array.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp" // IWYU pragma: keep

#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"


#include "gtest/gtest.h"
#include <core/numerics/ion_updater/ion_updater_def.hpp>
#include <cstddef>


namespace PHARE::core
{
// RUNTIME ENV VAR OVERRIDES
auto static const seeder    = get_env_as("PHARE_SEED", std::size_t{0}); // 0 == NO SEED!
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{10});
auto static const ppc       = get_env_as("PHARE_PPC", std::size_t{500});
auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{1});
auto static const dt        = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const do_cmp    = get_env_as("PHARE_COMPARE", std::size_t{1});
auto static const cmp_only  = get_env_as("PHARE_CMP_ONLY", std::size_t{0});
auto static const ref_only  = get_env_as("PHARE_REF_ONLY", std::size_t{0});

bool static const premain = []() { return true; }();


template<std::uint16_t impl, typename Ions>
auto get_updater_for(Ions const&)
{
    if constexpr (impl == 0)
        return IonUpdaterProxy<IonUpdater0<Ions>>();
    if constexpr (impl == 1)
        return IonUpdaterProxy<IonUpdater1<Ions>>();
    if constexpr (impl == 2)
        return IonUpdaterProxy<IonUpdater2<Ions>>();
}


template<typename Patches>
void ref_update(UpdaterMode mode, Patches& patches)
{
    using GridLayout_t = Patches::value_type::GridLayout_t;
    using Particles    = Patches::value_type::ParticleArray_t;
    using Boxing_t     = UpdaterSelectionBoxing<GridLayout_t, Particles>;

    for (auto& [layout, ions, _, electromag] : patches)
    {
        auto const domainBox = layout.AMRBox();
        auto const ghostBox  = grow(domainBox, GridLayout_t::nbrParticleGhosts());
        auto updater         = get_updater_for<0>(*ions);

        Boxing_t const boxing{layout, remove(ghostBox, domainBox)};
        updater.updatePopulations(mode, *ions, electromag, boxing, dt);
    }
}

template<std::uint16_t impl, typename Patches>
void cmp_update(UpdaterMode mode, Patches& patches)
{
    using GridLayout_t = Patches::value_type::GridLayout_t;
    using Particles    = Patches::value_type::ParticleArray_t;
    using Boxing_t     = UpdaterSelectionBoxing<GridLayout_t, Particles>;

    for (auto& [layout, ions, _, electromag] : patches)
    {
        auto const domainBox = layout.AMRBox();
        auto const ghostBox  = grow(domainBox, GridLayout_t::nbrParticleGhosts());
        auto updater         = get_updater_for<impl>(*ions);

        Boxing_t const boxing{layout, remove(ghostBox, domainBox)};
        updater.updatePopulations(mode, *ions, electromag, boxing, dt);
    }
}


template<typename Particles_t, typename GridLayout_t>
auto make_ions(GridLayout_t const& layout)
{
    auto constexpr static interp = GridLayout_t::interp_order;

    UsableIons<Particles_t, interp> ions{layout};
    std::optional<int> seed = std::nullopt;
    if (seeder > 0)
        seed = seeder;

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0ull);

    auto const disperse = [&](auto& particles) {
        delta_disperse(particles.domain_particles, seed);
        vary_velocity(particles.domain_particles, -6, 6, seed);

        delta_disperse(particles.level_ghost_particles, seed);
        vary_velocity(particles.level_ghost_particles, -6, 6, seed);
    };

    auto const particle_box = layout.AMRBox();

    auto add_particles = [&](auto& particles) {
        particles.domain_particles.reserve(particle_box.size() * ppc);
        add_particles_in(particles.domain_particles, particle_box, ppc);
    };

    add_particles(ions.populations[0].particles);
    disperse(ions.populations[0].particles);


    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), particle_box.size() * ppc);

    return ions;
}



template<typename Particles_t, typename GridLayout_t, typename Ions>
auto from_ions(GridLayout_t const& layout, Ions const& from)
{
    auto constexpr static interp = GridLayout_t::interp_order;

    UsableIons<Particles_t, interp> ions{layout, "protons"};
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0ull);

    auto _add_particles_from = [&](auto& src, auto& dst) { dst = src; };

    _add_particles_from.template operator()(from.populations[0].particles.domain_particles,
                                            ions.populations[0].particles.domain_particles);

    _add_particles_from.template operator()(from.populations[0].particles.level_ghost_particles,
                                            ions.populations[0].particles.level_ghost_particles);

    auto const particle_box = layout.AMRBox();
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), particle_box.size() * ppc);
    return std::move(ions);
}


template<typename GridLayout_t, typename P0, typename P1>
void check_particles(GridLayout_t const& layout, P0& ref, P1& cmp, double const atol)
{
    auto const box = layout.AMRBox();
    sort_particles(layout, ref);
    sort_particles(layout, cmp);

    EXPECT_EQ(ref.size(), cmp.size());

    auto const report = compare_particles(ref, cmp, atol);
    PHARE_LOG_LINE_STR("results: " << report.why());
    PHARE_LOG_LINE_STR("eg: " << ref[0]);
    EXPECT_TRUE(report);
}

template<typename GridLayout_t, typename R, typename C>
void compare(GridLayout_t const& layout, R& ref, C& cmp)
{
    double diff = 1e-14;

    check_particles(layout, ref.populations[0].particles.domain_particles,
                    cmp.populations[0].particles.domain_particles, diff);

    auto const freport = compare_tensor_fields(ref.populations[0].F, cmp.populations[0].F, diff);
    PHARE_LOG_LINE_STR("results: " << freport.why());
    EXPECT_TRUE(freport);

    auto const rhoCport = compare_fields(ref.populations[0].rhoC, cmp.populations[0].rhoC, diff);
    PHARE_LOG_LINE_STR("results: " << rhoCport.why());
    EXPECT_TRUE(rhoCport);

    auto const rhoPport = compare_fields(ref.populations[0].rhoP, cmp.populations[0].rhoP, diff);
    PHARE_LOG_LINE_STR("results: " << rhoPport.why());
    EXPECT_TRUE(rhoPport);
}

template<auto opts_, std::uint16_t impl_>
struct TestParam
{
    auto static constexpr opts         = opts_;
    auto static constexpr impl         = impl_;
    static constexpr auto dimension    = opts.dimension;
    static constexpr auto interp_order = opts.interp_order;
};




template<typename Param>
struct UpdatersComparisonTest : public ::testing::Test
{
    auto constexpr static opts   = Param::opts;
    auto constexpr static dim    = opts.dimension;
    auto constexpr static interp = opts.interp_order;

    using PHARE_Types        = core::PHARE_Types<opts>;
    using GridLayout_t       = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using RefParticleArray_t = ParticleArray<dim>;
    using CmpParticleArray_t = RefParticleArray_t;
    using UsableElectromag_t = UsableElectromag<dim>;
    using RefIons_t          = UsableIons<RefParticleArray_t, interp>;
    using RefEM_t            = UsableElectromag_t;
    using CmpIons_t          = UsableIons<CmpParticleArray_t, interp>;
    using CmpEM_t            = UsableElectromag_t;

    GridLayout_t const layout{cells};
    UpdaterMode const updater_mode = UpdaterMode::all;

    UpdatersComparisonTest() {}

    void run()
    {
        cmp_patches.reserve(n_patches);
        auto& ref = ref_patches.emplace_back(layout, make_ions<RefParticleArray_t>(layout));

        if (!ref_only)
            for (std::size_t i = 0; i < n_patches; i++)
                cmp_patches.emplace_back(layout, from_ions<CmpParticleArray_t>(layout, ref.ions));

        if (!cmp_only)
            ref_update(updater_mode, ref_patches);

        if (!ref_only)
            cmp_update<Param::impl>(updater_mode, cmp_patches);

        if (do_cmp)
            for (auto& cmp : cmp_patches)
                compare(*layout, ref_patches[0].ions, cmp.ions);
    }


    template<typename Ions_t, typename EM>
    struct Patch
    {
        using ParticleArray_t = Ions_t::particle_array_type;
        using GridLayout_t    = UpdatersComparisonTest<Param>::GridLayout_t;
        using Electromag_t    = EM::Super;

        GridLayout_t layout;
        Ions_t ions;
        EM em{layout};
        Electromag_t electromag = *em;

        std::string patchID() const { return "patch_id"; }
    };


    using RefPatch = Patch<RefIons_t, RefEM_t>;
    using CmpPatch = Patch<CmpIons_t, CmpEM_t>;

    std::vector<aggregate_adapter<RefPatch>> ref_patches{};
    std::vector<aggregate_adapter<CmpPatch>> cmp_patches{};
};

// clang-format off
using Permutations_t = ::testing::Types<
    TestParam<PHARE::SimOpts{1, 1}, 0>
  , TestParam<PHARE::SimOpts{2, 1}, 0>
  , TestParam<PHARE::SimOpts{3, 1}, 0>
  , TestParam<PHARE::SimOpts{1, 1}, 1>
  , TestParam<PHARE::SimOpts{2, 1}, 1>
  , TestParam<PHARE::SimOpts{3, 1}, 1>
  , TestParam<PHARE::SimOpts{1, 1}, 2>
  , TestParam<PHARE::SimOpts{2, 1}, 2>
  , TestParam<PHARE::SimOpts{3, 1}, 2>
>;
// clang-format on

TYPED_TEST_SUITE(UpdatersComparisonTest, Permutations_t, );

TYPED_TEST(UpdatersComparisonTest, updater_domain_only)
{
    this->run();
}

template<std::size_t dim> // used by gtest
void PrintTo(ParticleArray<dim> const& arr, std::ostream* os)
{
    *os << arr;
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
