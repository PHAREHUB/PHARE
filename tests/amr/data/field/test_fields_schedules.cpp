
#include "phare_core.hpp"
#include "phare/phare.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "amr/utilities/box/amr_box.hpp"

#include "tests/amr/test_hierarchy_fixtures.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include <core/data/ndarray/ndarray_vector.hpp>
#include <core/utilities/box/box.hpp>
#include <core/utilities/types.hpp>


#include "gtest/gtest.h"

#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/hier/HierarchyNeighbors.h>


static constexpr std::size_t ppc = 100;

template<std::size_t _dim, std::size_t interp>
struct TestParam
{
    auto constexpr static dim = _dim;
    using PhareTypes          = PHARE::core::PHARE_Types<dim, interp>;
    using GridLayout_t        = TestGridLayout<typename PhareTypes::GridLayout_t>;
    using Box_t               = PHARE::core::Box<int, dim>;
    using ParticleArray_t     = PHARE::core::ParticleArray<dim>;

    using Hierarchy_t
        = AfullHybridBasicHierarchy<dim, interp, defaultNbrRefinedParts<dim, interp>()>;
};

template<typename TestParam_>
struct FieldScheduleHierarchyTest : public ::testing::Test
{
    using TestParam           = TestParam_;
    using Hierarchy_t         = typename TestParam::Hierarchy_t;
    using ResourceManager_t   = typename Hierarchy_t::ResourcesManagerT;
    auto constexpr static dim = TestParam::dim;

    std::string const configFile
        = "test_fields_schedules_inputs/" + std::to_string(dim) + "d_L0.txt";
    Hierarchy_t hierarchy{configFile};
};

using FieldDatas = testing::Types</*TestParam<1, 1>, TestParam<1, 2>, TestParam<1, 3>,*/
                                  TestParam<2, 1> /*, TestParam<2, 2>, TestParam<2, 3>*/>;


TYPED_TEST_SUITE(FieldScheduleHierarchyTest, FieldDatas);
namespace PHARE::core
{



TYPED_TEST(FieldScheduleHierarchyTest, testing_hyhy_schedules)
{
    auto constexpr static dim    = TypeParam::dim;
    using GridLayout_t           = TestFixture::TestParam::GridLayout_t;
    using FieldData_t            = TestFixture::ResourceManager_t::UserField_t::patch_data_type;
    auto constexpr static interp = GridLayout_t::interp_order;
    auto constexpr static ghost_cells = GridLayout_t::nbrGhosts();

    auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& rm   = *this->hierarchy.resourcesManagerHybrid;
    auto& ions = this->hierarchy.hybridModel->state.ions;

    for (auto& patch : *lvl0)
    {
        auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);

        auto dataOnPatch = rm.setOnPatch(*patch, ions);
        for (auto& pop : ions)
        {
            auto const domGhostBox      = layout.AMRGhostBoxFor(pop.flux()[0].physicalQuantity());
            auto const primal_box       = shrink(domGhostBox, ghost_cells);
            auto const inner_primal_box = shrink(primal_box, 1);

            for (auto const& bix : layout.AMRToLocal(primal_box))
            {
                for (auto& f : pop.flux())
                    f(bix) = .2;
                pop.density()(bix) = .2;
            }
        }
    }

    this->hierarchy.messenger->fillFluxBorders(ions, *lvl0, 0);
    this->hierarchy.messenger->fillDensityBorders(ions, *lvl0, 0);


    auto proton_flux_x_id = *rm.getID("protons_flux_x");

    for (auto& patch : *lvl0)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);

        auto const field_data = SAMRAI_SHARED_PTR_CAST<FieldData_t, SAMRAI::hier::PatchData>(
            patch->getPatchData(proton_flux_x_id));

        SAMRAI::hier::HierarchyNeighbors const hier_nbrs{
            *this->hierarchy.basicHierarchy->hierarchy(), patch->getPatchLevelNumber(),
            patch->getPatchLevelNumber()};

        auto domainSamBox    = patch->getBox();
        auto const domainBox = phare_box_from<dim>(domainSamBox);

        auto const neighbors = core::generate(
            [](auto const& el) { return phare_box_from<dim>(el); },
            hier_nbrs.getSameLevelNeighbors(domainSamBox, patch->getPatchLevelNumber()));
        auto const ncells = core::sum_from(
            neighbors, [&](auto& el) { return (*(grow(el, 1) * domainBox)).size(); });

        if (mpi::rank() == 0)
            for (auto& pop : ions)
            {
                if (pop.name() == "protons")
                {
                    auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);

                    std::vector<double> vexpected(pop.flux()[0].size(), 1);
                    auto expected = core::make_array_view(vexpected.data(), pop.flux()[0].shape());

                    // auto const domGhostBox
                    //     = layout.AMRGhostBoxFor(pop.flux()[0].physicalQuantity());

                    // for (auto const& neighbor : neighbors)
                    //     if (auto p2poverlap = (domGhostBox * grow(neighbor, ghost_cells)))
                    //         for (auto const& bix : layout.AMRToLocal(*p2poverlap))
                    //             expected(*bix) += 1;

                    auto const domGhostBox
                        = layout.AMRGhostBoxFor(pop.flux()[0].physicalQuantity());
                    auto const primalDomainBox = [&]() {
                        auto box = domainBox;
                        box.upper += 1;
                        return box;
                    }();
                    auto const noverlap = shrink(primalDomainBox, 1 + (interp > 1));

                    for (auto const ghost_layer : domGhostBox.remove(noverlap))
                        for (auto const& neighbor : neighbors)
                        {
                            auto nghostbox = grow(neighbor, GridLayout_t::nbrGhosts());
                            nghostbox.upper += 1;
                            if (auto p2poverlap = (nghostbox * ghost_layer))
                                for (auto const& bix : layout.AMRToLocal(*p2poverlap))
                                    expected(*bix) += 1;
                        }

                    if (vexpected != field_data->field.vector())
                        field_data->gridLayout.evalOnGhostBox(pop.flux()[0], [&](auto... args) {
                            PHARE_LOG_LINE_SS((Point{args...}.as_signed() + domGhostBox.lower)
                                              << " " << pop.flux()[0](args...));
                        });

                    for (auto const& e : field_data->field)
                        EXPECT_NE(e, 0); // :/
                }
            }
    }
}



} // namespace PHARE::core


int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
