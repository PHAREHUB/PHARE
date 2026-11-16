// The fully kinetic J ghost fill copies every ghost node from the patch that owns it, corners
// and edges included: the current filter reads them all.

#include "phare_mpi.hpp"

#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "tests/simulator/per_test.hpp"

#include <cmath>
#include <limits>
#include <numbers>


using namespace PHARE;

static constexpr SimOpts opts{
    .dimension = 2, .interp_order = 1, .nbRefinedPart = 4, .fully_kinetic = true};
using Simulator_t = SimulatorTestParam<opts>;
using FKModel     = solver::PHARE_Types<opts>::FK::Model_t;
using GridLayout  = FKModel::gridlayout_type;


TEST(FullyKineticCurrentGhosts, areCopiedFromTheirOwnerIncludingCorners)
{
    Simulator_t sim;
    auto& model     = *sim.getFKModel();
    auto& rm        = *model.resourcesManager;
    auto& J         = model.state.J;
    auto& level     = *sim.hierarchy->getPatchLevel(0);
    auto& messenger = dynamic_cast<amr::HybridMessenger<FKModel>&>(
        sim.getMultiPhysicsIntegrator()->messenger(0));

    double constexpr L = 16 * 0.04;
    double constexpr k = 2 * std::numbers::pi / L;
    auto const nan     = std::numeric_limits<double>::quiet_NaN();

    // periodic, different on each component
    auto const f = [&](std::size_t c, auto const& r) {
        return std::sin(k * r[0] + c) * std::cos(k * r[1] - 0.3 * c) + c;
    };

    auto const isPhysical = [](auto const& layout, auto const& Jc, auto const&... ijk) {
        std::array<std::uint32_t, 2> const index{static_cast<std::uint32_t>(ijk)...};
        for (std::size_t d = 0; d < 2; ++d)
        {
            auto const direction = static_cast<core::Direction>(d);
            if (index[d] < layout.physicalStartIndex(Jc, direction)
                or index[d] > layout.physicalEndIndex(Jc, direction))
                return false;
        }
        return true;
    };

    auto const position = [](auto const& layout, auto const& Jc, auto const&... ijk) {
        return layout.fieldNodeCoordinates(
            Jc, layout.localToAMR(core::Point{static_cast<std::uint32_t>(ijk)...}.as_signed()));
    };

    std::size_t nbrPatches = 0;
    for (auto& patch : rm.enumerate(level, J))
    {
        ++nbrPatches;
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        for (std::size_t c = 0; c < 3; ++c)
        {
            auto& Jc = J[c];
            layout.evalOnGhostBox(Jc, [&](auto const&... ijk) {
                Jc(ijk...)
                    = isPhysical(layout, Jc, ijk...) ? f(c, position(layout, Jc, ijk...)) : nan;
            });
        }
    }
    ASSERT_EQ(nbrPatches, 4u);

    messenger.fillCurrentGhosts(J, level, 0.);

    std::size_t corners = 0;
    for (auto& patch : rm.enumerate(level, J))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        for (std::size_t c = 0; c < 3; ++c)
        {
            auto& Jc = J[c];
            layout.evalOnGhostBox(Jc, [&](auto const&... ijk) {
                auto const expected = f(c, position(layout, Jc, ijk...));
                EXPECT_NEAR(Jc(ijk...), expected, 1e-14) << "component " << c;

                std::array<std::uint32_t, 2> const index{static_cast<std::uint32_t>(ijk)...};
                std::size_t outside = 0;
                for (std::size_t d = 0; d < 2; ++d)
                {
                    auto const direction = static_cast<core::Direction>(d);
                    outside += index[d] < layout.physicalStartIndex(Jc, direction)
                               or index[d] > layout.physicalEndIndex(Jc, direction);
                }
                corners += outside == 2;
            });
        }
    }
    EXPECT_GT(corners, 0u);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
