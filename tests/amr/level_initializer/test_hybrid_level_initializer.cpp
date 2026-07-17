

#include "phare_core.hpp"
#include "phare_solver.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/field/field_box.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "amr/level_initializer/hybrid_level_initializer.hpp"

#include "tests/amr/amr.hpp"
#include "tests/amr/test_hierarchy_fixtures.hpp"
#include "tests/core/data/field/test_field_fixtures.hpp"
#include "tests/core/data/grid/test_grid_fixtures.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

namespace PHARE::amr
{

// Minimal tag strategy: only allocates and registers — no model.initialize.
template<typename HybridModel>
class AllocOnlyTagStrategy : public SAMRAI::mesh::StandardTagAndInitStrategy
{
    using HybridMessenger_t = HybridMessenger<HybridModel>;

public:
    AllocOnlyTagStrategy(std::shared_ptr<HybridModel> model,
                         std::shared_ptr<SolverPPC<HybridModel, SAMRAI_Types>> solver,
                         std::shared_ptr<HybridMessenger_t> messenger)
        : model_{std::move(model)}
        , solver_{std::move(solver)}
        , messenger_{std::move(messenger)}
    {
        auto infoFromFiner   = messenger_->emptyInfoFromFiner();
        auto infoFromCoarser = messenger_->emptyInfoFromCoarser();
        model_->fillMessengerInfo(infoFromFiner);
        model_->fillMessengerInfo(infoFromCoarser);
        solver_->fillMessengerInfo(infoFromFiner);
        messenger_->registerQuantities(std::move(infoFromCoarser), std::move(infoFromFiner));
    }

    void initializeLevelData(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                             int const levelNumber, double const initDataTime, bool const,
                             bool const,
                             std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*oldLevel*/
                             = std::shared_ptr<SAMRAI::hier::PatchLevel>(),
                             bool const allocateData = true) override
    {
        auto level = hierarchy->getPatchLevel(levelNumber);
        if (allocateData)
            for (auto patch : *level)
            {
                model_->allocate(*patch, initDataTime);
                solver_->allocate(*model_, *patch, initDataTime);
                messenger_->allocate(*patch, initDataTime);
            }
        messenger_->registerLevel(hierarchy, levelNumber);
    }

    void resetHierarchyConfiguration(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const&,
                                     int const, int const) override
    {
    }

private:
    std::shared_ptr<HybridModel> model_;
    std::shared_ptr<SolverPPC<HybridModel, SAMRAI_Types>> solver_;
    std::shared_ptr<HybridMessenger_t> messenger_;
};


template<SimOpts opts>
struct LevelInitRunner
{
    auto constexpr static dim = opts.dimension;

    using Solver_t     = solver::PHARE_Types<opts>;
    using GridLayout_t = Solver_t::GridLayout_t;
    using HybridModelT = Solver_t::HybridModel_t;

    using HybridHybridT
        = HybridHybridMessengerStrategy<HybridModelT, typename Solver_t::RefinementParams>;

    std::string const configFile
        = "test_hybrid_level_initializer_inputs/" + std::to_string(dim) + "d_config.txt";

    PHARE::initializer::PHAREDict dict{createDict<dim>()};

    std::shared_ptr<typename HybridModelT::resources_manager_type> rm{
        std::make_shared<typename HybridModelT::resources_manager_type>()};
    std::shared_ptr<HybridModelT> hybridModel{std::make_shared<HybridModelT>(dict, rm)};

    std::shared_ptr<HybridMessenger<HybridModelT>> messenger{
        std::make_shared<HybridMessenger<HybridModelT>>(std::make_unique<HybridHybridT>(rm, 0))};

    std::shared_ptr<SolverPPC<HybridModelT, SAMRAI_Types>> solver{
        std::make_shared<SolverPPC<HybridModelT, SAMRAI_Types>>(dict["simulation"]["algo"])};

    bool const _registered = [&] {
        hybridModel->resourcesManager->registerResources(*hybridModel);
        solver->registerResources(*hybridModel);
        return true;
    }();

    std::shared_ptr<AllocOnlyTagStrategy<HybridModelT>> tagStrat{
        std::make_shared<AllocOnlyTagStrategy<HybridModelT>>(hybridModel, solver, messenger)};

    std::shared_ptr<TestIntegratorStrat> integratorStrat{std::make_shared<TestIntegratorStrat>()};

    std::shared_ptr<BasicHierarchy> basicHierarchy{
        std::make_shared<BasicHierarchy>(2, dim, tagStrat.get(), integratorStrat, configFile)};

    solver::HybridLevelInitializer<HybridModelT> levelInit{core::OhmInfo{0.0, 0.0001}};

    void initialize()
    {
        auto lvl0 = basicHierarchy->hierarchy()->getPatchLevel(0);
        for (auto& patch : rm->enumerate(*lvl0, *hybridModel))
        {
            auto const layout = layoutFromPatch<GridLayout_t>(*patch);
            for (auto& pop : hybridModel->state.ions)
                for (auto& c : hybridModel->state.electromag.B)
                    c.fill(1);
        }
        levelInit.initialize(basicHierarchy->hierarchy(), 0, nullptr, *hybridModel, *messenger, 0.0,
                             false);
    }

    // Per-patch snapshot of named scalar fields as plain doubles.
    // Survives runner destruction so two runners can be compared sequentially.
    struct PatchSnapshot
    {
        core::Box<int, dim> amr_box;
        std::map<std::string, std::vector<double>> fields; // name -> domain-cell values
    };

    void add_field(PatchSnapshot& ps, std::string const& name, auto const& field,
                   GridLayout_t const& layout)
    {
        auto const reduced = core::reduce_single(field);
        auto& values       = ps.fields[name];
        for (auto const& idx : layout.domainBoxFor(reduced.physicalQuantity()))
            values.push_back(static_cast<double>(reduced(idx)));
    }

    void add_vecfield(PatchSnapshot& ps, std::string const& name, auto const& vf,
                      GridLayout_t const& layout)
    {
        static constexpr std::array comps{"x", "y", "z"};
        for (int c = 0; c < 3; ++c)
            add_field(ps, name + "_" + comps[c], vf[c], layout);
    }

    // Collect B, E, J, ion density, ion velocity, and per-pop density+flux on lvl0.
    auto collect_fields_lvl0()
    {
        std::vector<PatchSnapshot> result;
        auto lvl0   = basicHierarchy->hierarchy()->getPatchLevel(0);
        auto& state = hybridModel->state;

        for (auto& patch : rm->enumerate(*lvl0, *hybridModel))
        {
            auto const layout = layoutFromPatch<GridLayout_t>(*patch);
            PatchSnapshot ps;
            ps.amr_box = layout.AMRBox();

            add_vecfield(ps, "B", state.electromag.B, layout);
            add_vecfield(ps, "E", state.electromag.E, layout);
            add_vecfield(ps, "J", state.J, layout);
            add_field(ps, "Ni", state.ions.chargeDensity(), layout);
            add_vecfield(ps, "Vi", state.ions.velocity(), layout);

            for (auto const& pop : state.ions)
            {
                add_field(ps, pop.name() + "_density", pop.particleDensity(), layout);
                add_vecfield(ps, pop.name() + "_flux", pop.flux(), layout);
            }

            result.push_back(std::move(ps));
        }

        std::sort(result.begin(), result.end(), [](PatchSnapshot const& a, PatchSnapshot const& b) {
            for (int d = 0; d < dim; ++d)
                if (a.amr_box.lower[d] != b.amr_box.lower[d])
                    return a.amr_box.lower[d] < b.amr_box.lower[d];
            return false;
        });
        return result;
    }
};


// ---- typed tests: per-layout sanity checks --------------------------------

template<SimOpts opts>
struct TestParam
{
    auto constexpr static dim = opts.dimension;
    using Runner_t            = LevelInitRunner<opts>;
};

template<typename TestParam_>
struct HybridLevelInitTest : public ::testing::Test
{
    using Runner_t = TestParam_::Runner_t;
    Runner_t runner;
};

// clang-format off
using TestTypes = testing::Types<
    TestParam<SimOpts{.dimension=2}>
   ,TestParam<SimOpts{.dimension=3}>
PHARE_WITH_MKN_GPU(
   ,TestParam<SimOpts{.dimension=2, .layout_mode=LayoutMode::AoSTS}>
   ,TestParam<SimOpts{.dimension=3, .layout_mode=LayoutMode::AoSTS}>
)
>;
// clang-format on

TYPED_TEST_SUITE(HybridLevelInitTest, TestTypes, );


TYPED_TEST(HybridLevelInitTest, E_domain_no_nan_after_init)
{
    using GridLayout_t = TestFixture::Runner_t::GridLayout_t;
    this->runner.initialize();

    auto lvl0 = this->runner.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& em  = this->runner.hybridModel->state.electromag;

    for (auto& patch : this->runner.rm->enumerate(*lvl0, *this->runner.hybridModel))
    {
        auto const layout = layoutFromPatch<GridLayout_t>(*patch);
        for (auto& field : em.E)
        {
            auto const reduced = core::reduce_single(field);
            for (auto const& idx : layout.domainBoxFor(reduced.physicalQuantity()))
                EXPECT_FALSE(std::isnan(reduced(idx)))
                    << "NaN in E." << field.name() << " at " << idx << " patch=" << layout.AMRBox()
                    << " rank=" << mpi::rank();
        }
    }
}


TYPED_TEST(HybridLevelInitTest, E_init_per_patch_log)
{
    using GridLayout_t = TestFixture::Runner_t::GridLayout_t;
    this->runner.initialize();

    auto lvl0 = this->runner.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& em  = this->runner.hybridModel->state.electromag;

    for (auto& patch : this->runner.rm->enumerate(*lvl0, *this->runner.hybridModel))
    {
        auto const layout = layoutFromPatch<GridLayout_t>(*patch);
        PHARE_LOG_LINE_SS("rank=" << mpi::rank() << " patch=" << layout.AMRBox());
        for (auto& field : em.E)
        {
            auto const reduced = core::reduce_single(field);
            PHARE_LOG_LINE_SS("  " << field.name() << " sum=" << core::sum_field(reduced)
                                   << " sum_not_nan=" << core::sum_not_nan(reduced));
        }
    }
}


TYPED_TEST(HybridLevelInitTest, density_and_flux_nonzero_after_init)
{
    using GridLayout_t = TestFixture::Runner_t::GridLayout_t;
    this->runner.initialize();

    auto lvl0  = this->runner.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& ions = this->runner.hybridModel->state.ions;

    for (auto& patch : this->runner.rm->enumerate(*lvl0, *this->runner.hybridModel))
    {
        auto const layout = layoutFromPatch<GridLayout_t>(*patch);

        for (auto const& pop : ions)
        {
            auto const rho = core::reduce_single(pop.particleDensity());
            PHARE_LOG_LINE_SS("rank=" << mpi::rank() << " patch=" << layout.AMRBox() << " "
                                      << pop.name()
                                      << "_density sum_not_nan=" << core::sum_not_nan(rho));
            for (auto const& idx : layout.domainBoxFor(rho.physicalQuantity()))
                EXPECT_GT(std::abs(rho(idx)), .01)
                    << pop.name() << "_density near-zero at " << idx << " val=" << rho(idx)
                    << " patch=" << layout.AMRBox() << " rank=" << mpi::rank();

            for (int c = 0; c < 3; ++c)
            {
                auto const flux = core::reduce_single(pop.flux()[c]);
                PHARE_LOG_LINE_SS("rank=" << mpi::rank() << " patch=" << layout.AMRBox() << " "
                                          << pop.name() << "_flux[" << c
                                          << "] sum_not_nan=" << core::sum_not_nan(flux));
                for (auto const& idx : layout.domainBoxFor(flux.physicalQuantity()))
                    EXPECT_GT(std::abs(flux(idx)), .01)
                        << pop.name() << "_flux[" << c << "] near-zero at " << idx
                        << " val=" << flux(idx) << " patch=" << layout.AMRBox()
                        << " rank=" << mpi::rank();
            }
        }
    }
}


// // ---- comparison: AoSMapped (reference) vs AoSTS ---------------------------

#ifdef PHARE_WITH_MKN_GPU

TEST(HybridLevelInitComparison, fields_AoSMapped_vs_AoSTS_2d)
{
    constexpr double atol = 1e-12;
    using RefRunner       = LevelInitRunner<SimOpts{.dimension = 2}>;
    using CmpRunner = LevelInitRunner<SimOpts{.dimension = 2, .layout_mode = LayoutMode::AoSTS}>;

    std::vector<RefRunner::PatchSnapshot> ref_snaps;
    {
        RefRunner ref;
        ref.initialize();
        ref_snaps = ref.collect_fields_lvl0();
    }

    std::vector<CmpRunner::PatchSnapshot> cmp_snaps;
    {
        CmpRunner cmp;
        cmp.initialize();
        cmp_snaps = cmp.collect_fields_lvl0();
    }

    ASSERT_EQ(ref_snaps.size(), cmp_snaps.size()) << "patch count mismatch rank=" << mpi::rank();

    for (std::size_t pi = 0; pi < ref_snaps.size(); ++pi)
    {
        auto const& rsnap = ref_snaps[pi];
        auto const& csnap = cmp_snaps[pi];

        ASSERT_EQ(rsnap.amr_box, csnap.amr_box)
            << "patch box mismatch at index " << pi << " rank=" << mpi::rank();

        for (auto const& [name, rv] : rsnap.fields)
        {
            auto it = csnap.fields.find(name);
            ASSERT_NE(it, csnap.fields.end()) << "field " << name << " missing in cmp";
            auto const& cv = it->second;
            ASSERT_EQ(rv.size(), cv.size())
                << "size mismatch for " << name << " patch=" << rsnap.amr_box;
            for (std::size_t i = 0; i < rv.size(); ++i)
                EXPECT_NEAR(rv[i], cv[i], atol)
                    << name << "[" << i << "] ref=" << rv[i] << " cmp=" << cv[i]
                    << " patch=" << rsnap.amr_box << " rank=" << mpi::rank();
        }
    }
}

TEST(HybridLevelInitComparison, fields_AoSMapped_vs_AoSTS_3d)
{
    constexpr double atol = 1e-12;
    using RefRunner       = LevelInitRunner<SimOpts{.dimension = 3}>;
    using CmpRunner = LevelInitRunner<SimOpts{.dimension = 3, .layout_mode = LayoutMode::AoSTS}>;

    std::vector<RefRunner::PatchSnapshot> ref_snaps;
    {
        RefRunner ref;
        ref.initialize();
        ref_snaps = ref.collect_fields_lvl0();
    }

    std::vector<CmpRunner::PatchSnapshot> cmp_snaps;
    {
        CmpRunner cmp;
        cmp.initialize();
        cmp_snaps = cmp.collect_fields_lvl0();
    }

    ASSERT_EQ(ref_snaps.size(), cmp_snaps.size()) << "patch count mismatch rank=" << mpi::rank();

    for (std::size_t pi = 0; pi < ref_snaps.size(); ++pi)
    {
        auto const& rsnap = ref_snaps[pi];
        auto const& csnap = cmp_snaps[pi];

        ASSERT_EQ(rsnap.amr_box, csnap.amr_box)
            << "patch box mismatch at index " << pi << " rank=" << mpi::rank();

        for (auto const& [name, rv] : rsnap.fields)
        {
            auto it = csnap.fields.find(name);
            ASSERT_NE(it, csnap.fields.end()) << "field " << name << " missing in cmp";
            auto const& cv = it->second;
            ASSERT_EQ(rv.size(), cv.size())
                << "size mismatch for " << name << " patch=" << rsnap.amr_box;
            for (std::size_t i = 0; i < rv.size(); ++i)
                EXPECT_NEAR(rv[i], cv[i], atol)
                    << name << "[" << i << "] ref=" << rv[i] << " cmp=" << cv[i]
                    << " patch=" << rsnap.amr_box << " rank=" << mpi::rank();
        }
    }
}

#endif // PHARE_WITH_MKN_GPU


} // namespace PHARE::amr


int main(int argc, char** argv)
{
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
