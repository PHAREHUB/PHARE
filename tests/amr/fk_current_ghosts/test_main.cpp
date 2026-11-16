// The fully kinetic J and rho ghost fills copy every ghost node from the patch that owns it,
// corners and edges included: the filter reads them all.

#include "phare_mpi.hpp"

#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "core/data/field/field_traits.hpp"

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


// the components of a vector field, or the field itself
template<typename Resource>
auto components(Resource& resource)
{
    if constexpr (core::IsField<Resource>)
        return std::array{&resource};
    else
        return std::array{&resource[0], &resource[1], &resource[2]};
}


// sets the physical nodes of each component to a periodic function and the ghosts to NaN,
// fills the ghosts, then every ghost node must hold the function, corners included
template<typename Resource, typename FillGhosts>
void checkGhostsAreCopiedFromTheirOwner(Simulator_t& sim, Resource& resource,
                                        FillGhosts&& fillGhosts)
{
    auto& model = *sim.getFKModel();
    auto& rm    = *model.resourcesManager;
    auto& level = *sim.hierarchy->getPatchLevel(0);

    double constexpr L = 16 * 0.04;
    double constexpr k = 2 * std::numbers::pi / L;
    auto const nan     = std::numeric_limits<double>::quiet_NaN();

    // periodic, different on each component
    auto const f = [&](std::size_t c, auto const& r) {
        return std::sin(k * r[0] + c) * std::cos(k * r[1] - 0.3 * c) + c;
    };

    auto const isPhysical = [](auto const& layout, auto const& field, auto const&... ijk) {
        std::array<std::uint32_t, 2> const index{static_cast<std::uint32_t>(ijk)...};
        for (std::size_t d = 0; d < 2; ++d)
        {
            auto const direction = static_cast<core::Direction>(d);
            if (index[d] < layout.physicalStartIndex(field, direction)
                or index[d] > layout.physicalEndIndex(field, direction))
                return false;
        }
        return true;
    };

    auto const position = [](auto const& layout, auto const& field, auto const&... ijk) {
        return layout.fieldNodeCoordinates(
            field, layout.localToAMR(core::Point{static_cast<std::uint32_t>(ijk)...}.as_signed()));
    };

    std::size_t nbrPatches = 0;
    for (auto& patch : rm.enumerate(level, resource))
    {
        ++nbrPatches;
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        std::size_t c     = 0;
        for (auto* field : components(resource))
        {
            layout.evalOnGhostBox(*field, [&](auto const&... ijk) {
                (*field)(ijk...) = isPhysical(layout, *field, ijk...)
                                       ? f(c, position(layout, *field, ijk...))
                                       : nan;
            });
            ++c;
        }
    }
    ASSERT_EQ(nbrPatches, 4u);

    fillGhosts(level);

    std::size_t corners = 0;
    for (auto& patch : rm.enumerate(level, resource))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        std::size_t c     = 0;
        for (auto* field : components(resource))
        {
            layout.evalOnGhostBox(*field, [&](auto const&... ijk) {
                auto const expected = f(c, position(layout, *field, ijk...));
                EXPECT_NEAR((*field)(ijk...), expected, 1e-14) << field->name();

                std::array<std::uint32_t, 2> const index{static_cast<std::uint32_t>(ijk)...};
                std::size_t outside = 0;
                for (std::size_t d = 0; d < 2; ++d)
                {
                    auto const direction = static_cast<core::Direction>(d);
                    outside += index[d] < layout.physicalStartIndex(*field, direction)
                               or index[d] > layout.physicalEndIndex(*field, direction);
                }
                corners += outside == 2;
            });
            ++c;
        }
    }
    EXPECT_GT(corners, 0u);
}


auto& messengerOf(Simulator_t& sim)
{
    return dynamic_cast<amr::HybridMessenger<FKModel>&>(
        sim.getMultiPhysicsIntegrator()->messenger(0));
}


TEST(FullyKineticGhosts, currentAndChargeDensityGhostsAreCopiedFromTheirOwnerWithCorners)
{
    Simulator_t sim;
    auto& state     = sim.getFKModel()->state;
    auto& messenger = messengerOf(sim);

    checkGhostsAreCopiedFromTheirOwner(
        sim, state.J, [&](auto& level) { messenger.fillCurrentGhosts(state.J, level, 0.); });
    checkGhostsAreCopiedFromTheirOwner(sim, state.rho, [&](auto& level) {
        messenger.fillChargeDensityGhosts(state.rho, level, 0.);
    });
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
