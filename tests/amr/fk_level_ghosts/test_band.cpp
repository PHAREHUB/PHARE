// The fully kinetic level ghost band must be wide enough (P4 plan, band width decision 11-17):
// the current filter reads J and rho on the first ghost layer outside a fine level, so every
// particle whose current or charge reaches that layer must be in the band. Here J and rho
// deposited on the fine level by depositFKCurrent are checked, on the physical nodes and on the
// first ghost layer outside the level, against a sum over all the particles of the plasma: the
// fine domain particles inside the level and, outside, the children of every coarse domain
// particle, split here with the splitter of the band. Each particle is in it once.

#include "phare_mpi.hpp"

#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/fk_current_deposit.hpp"
#include "amr/data/particles/refine/particles_data_split.hpp"
#include "core/data/field/field_traits.hpp"

#include "tests/simulator/per_test.hpp"

#include <set>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>


using namespace PHARE;

static constexpr std::size_t interp = PHARE_TEST_INTERP_ORDER;
static constexpr SimOpts opts{
    .dimension = 2, .interp_order = interp, .nbRefinedPart = 4, .fully_kinetic = true};
using Simulator_t = SimulatorTestParam<opts>;
using FKTypes     = solver::PHARE_Types<opts>::FK;
using FKModel     = FKTypes::Model_t;
using GridLayout  = FKModel::gridlayout_type;
using Key         = std::array<int, 2>;


// the hybrid widths as instantiated
template<std::size_t order, bool fk>
using Options
    = HybridFieldOptions<SimOpts{.dimension = 2, .interp_order = order, .fully_kinetic = fk}>;
static_assert(Options<1, false>::particle_ghost_width == 1
              and Options<1, false>::field_ghost_width == 2);
static_assert(Options<2, false>::particle_ghost_width == 2
              and Options<2, false>::field_ghost_width == 4);
static_assert(Options<3, false>::particle_ghost_width == 2
              and Options<3, false>::field_ghost_width == 4);


static constexpr double dtFine = 0.002 / 2; // job time step over the refinement ratio

// one simulator per process, destroyed before SAMRAI is shut down (main)
std::unique_ptr<Simulator_t> sim_;
Simulator_t& simulator()
{
    if (!sim_)
        sim_ = std::make_unique<Simulator_t>("job_o" + std::to_string(interp));
    return *sim_;
}


// the B-spline of the interpolation order, centered on the node, written here independently of
// the code under test
double shape(double x)
{
    x = std::abs(x);
    if constexpr (interp == 1)
        return x < 1 ? 1 - x : 0;
    else if constexpr (interp == 2)
        return x < 0.5 ? 0.75 - x * x : x < 1.5 ? 0.5 * (1.5 - x) * (1.5 - x) : 0;
    else
        return x < 1   ? 2. / 3. - x * x + 0.5 * x * x * x
               : x < 2 ? (2 - x) * (2 - x) * (2 - x) / 6.
                       : 0;
}


struct Charge
{
    std::array<int, 2> iCell; // fine AMR cell
    std::array<double, 2> delta;
    std::array<double, 3> v;
    double qw;
};


// local index of the node in [physical start - layers, physical end + layers] in both directions
bool within(GridLayout const& layout, auto const& field, std::uint32_t const layers, auto i, auto j)
{
    std::array<std::uint32_t, 2> const index{static_cast<std::uint32_t>(i),
                                             static_cast<std::uint32_t>(j)};
    for (std::size_t d = 0; d < 2; ++d)
    {
        auto const direction = static_cast<core::Direction>(d);
        if (index[d] + layers < layout.physicalStartIndex(field, direction)
            or index[d] > layout.physicalEndIndex(field, direction) + layers)
            return false;
    }
    return true;
}

Key amrIndex(GridLayout const& layout, auto i, auto j)
{
    auto const amr = layout.localToAMR(
        core::Point{static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j)}.as_signed());
    return {amr[0], amr[1]};
}


// J and rho of the particles on the local arrays of a patch: rho at x - dt/2 v, J by the
// Esirkepov formulas over x - dt/2 v -> x + dt/2 v (Esirkepov 2001, eqs. 31-38: the current along
// d is the running sum of the shape change along d times the mean shape across)
struct Reference
{
    std::array<std::size_t, 2> shape_; // of rho, primal in both directions
    std::vector<double> rho, Jx, Jy, Jz;

    explicit Reference(std::array<std::size_t, 2> shape)
        : shape_{shape}
        , rho(shape[0] * shape[1])
        , Jx(rho.size())
        , Jy(rho.size())
        , Jz(rho.size())
    {
    }

    std::size_t at(int i, int j) const { return i * shape_[1] + j; }
    bool inside(int i, int j) const
    {
        return i >= 0 and j >= 0 and i < static_cast<int>(shape_[0])
               and j < static_cast<int>(shape_[1]);
    }

    void add(GridLayout const& layout, Charge const& p)
    {
        auto const& dl               = layout.meshSize();
        double constexpr halfSupport = 0.5 * (interp + 1);
        int constexpr nodes          = interp + 4; // both shapes, a node of margin each side
        std::array<int, 2> lo        = {};
        std::array<std::array<double, nodes>, 2> S0{}, dS{};
        for (std::size_t d = 0; d < 2; ++d)
        {
            auto const pos = p.iCell[d] - layout.AMRBox().lower[d]
                             + static_cast<int>(layout.physicalStartIndex(
                                 core::QtyCentering::primal, static_cast<core::Direction>(d)))
                             + p.delta[d];
            auto const shift = 0.5 * dtFine * p.v[d] / dl[d];
            lo[d]            = static_cast<int>(std::ceil(pos - shift - halfSupport)) - 1;
            for (int n = 0; n < nodes; ++n)
            {
                S0[d][n] = shape(pos - shift - (lo[d] + n));
                dS[d][n] = shape(pos + shift - (lo[d] + n)) - S0[d][n];
            }
        }

        double const cx = -p.qw * dl[0] / dtFine, cy = -p.qw * dl[1] / dtFine;
        for (int a = 0; a < nodes; ++a)
            for (int b = 0; b < nodes; ++b)
            {
                int const i = lo[0] + a, j = lo[1] + b;
                if (!inside(i, j))
                    continue;
                rho[at(i, j)] += p.qw * S0[0][a] * S0[1][b];
                Jz[at(i, j)]
                    += p.qw * p.v[2]
                       * (S0[0][a] * S0[1][b] + 0.5 * (S0[0][a] * dS[1][b] + dS[0][a] * S0[1][b])
                          + dS[0][a] * dS[1][b] / 3.);
                double sumX = 0, sumY = 0; // dual node i + 1/2 (j + 1/2) has local index i (j)
                for (int k = 0; k <= a; ++k)
                    sumX += dS[0][k];
                for (int k = 0; k <= b; ++k)
                    sumY += dS[1][k];
                Jx[at(i, j)] += cx * sumX * (S0[1][b] + 0.5 * dS[1][b]);
                Jy[at(i, j)] += cy * sumY * (S0[0][a] + 0.5 * dS[0][a]);
            }
    }
};


TEST(FullyKineticLevelGhostBand, firstGhostLayerOutsideTheLevelGetsEveryParticle)
{
    auto& sim    = simulator();
    auto& model  = *sim.getFKModel();
    auto& rm     = *model.resourcesManager;
    auto& coarse = *sim.hierarchy->getPatchLevel(0);
    auto& fine   = *sim.hierarchy->getPatchLevel(1);
    auto& ions   = model.state.ions;
    auto& J      = model.state.J;
    auto& rho    = model.state.rho;
    ASSERT_EQ(sim.hierarchy->getNumberOfLevels(), 2);
    ASSERT_GT(model.divergence_cleaning, 0.);
    ASSERT_EQ(model.current_filter_passes, 0u);

    auto& messenger = dynamic_cast<amr::HybridMessenger<FKModel>&>(
        sim.getMultiPhysicsIntegrator()->messenger(1));
    solver::depositFKCurrent(fine, model, messenger, dtFine, 1., 1., 0.);

    auto const charge = [](auto const& pop, auto const& p) {
        return Charge{p.iCell,
                      {static_cast<double>(p.delta[0]), static_cast<double>(p.delta[1])},
                      p.v,
                      pop.charge() * p.weight};
    };

    // the fine level cells, and the nodes some fine patch owns (physical)
    std::vector<core::Box<int, 2>> fineBoxes;
    std::array<std::set<Key>, 4> owned; // Jx, Jy, Jz, rho
    for (auto& patch : rm.enumerate(fine, J, rho))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        fineBoxes.push_back(layout.AMRBox());
        std::size_t c = 0;
        for (auto* field : {&J[0], &J[1], &J[2], &rho})
        {
            layout.evalOnGhostBox(*field, [&](auto i, auto j) {
                if (within(layout, *field, 0, i, j))
                    owned[c].insert(amrIndex(layout, i, j));
            });
            ++c;
        }
    }
    auto const inLevel = [&](auto const& iCell) {
        return std::any_of(fineBoxes.begin(), fineBoxes.end(),
                           [&](auto const& box) { return core::isIn(core::Point{iCell}, box); });
    };

    std::vector<Charge> plasma;
    for (auto& patch : rm.enumerate(fine, ions))
        for (auto const& pop : ions)
            for (auto const& p : pop.domainParticles())
                plasma.push_back(charge(pop, p));
    auto const nbrFine = plasma.size();

    FKTypes::Splitter_t split;
    for (auto& patch : rm.enumerate(coarse, ions))
        for (auto const& pop : ions)
            for (auto const& p : pop.domainParticles())
            {
                std::array<std::decay_t<decltype(p)>, opts.nbRefinedPart> children;
                split(amr::toFineGrid<interp>(p), children);
                for (auto const& child : children)
                    if (!inLevel(child.iCell))
                        plasma.push_back(charge(pop, child));
            }
    ASSERT_GT(plasma.size(), nbrFine);

    std::array<double, 4> maxValue{}, maxError{};
    std::array<std::size_t, 4> nbrFirstLayer{};
    for (auto& patch : rm.enumerate(fine, J, rho))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        Reference reference{{rho.shape()[0], rho.shape()[1]}};
        for (auto const& p : plasma)
            reference.add(layout, p);

        std::size_t c = 0;
        for (auto const& [field, expected] :
             {std::pair{&J[0], &reference.Jx}, std::pair{&J[1], &reference.Jy},
              std::pair{&J[2], &reference.Jz}, std::pair{&rho, &reference.rho}})
        {
            layout.evalOnGhostBox(*field, [&](auto i, auto j) {
                if (!within(layout, *field, 1, i, j))
                    return;
                bool const physical = within(layout, *field, 0, i, j);
                if (!physical and owned[c].contains(amrIndex(layout, i, j)))
                    return; // first ghost layer inside the level, another patch's node
                auto const value = (*expected)[reference.at(i, j)];
                maxValue[c]      = std::max(maxValue[c], std::abs(value));
                maxError[c]      = std::max(maxError[c], std::abs((*field)(i, j) - value));
                nbrFirstLayer[c] += !physical;
            });
            ++c;
        }
    }

    std::size_t c = 0;
    for (auto const name : {"Jx", "Jy", "Jz", "rho"})
    {
        std::cout << "order " << interp << " " << name << ": " << nbrFirstLayer[c]
                  << " first ghost layer nodes outside the level, max error "
                  << maxError[c] / maxValue[c] << " of max |" << name << "|\n";
        EXPECT_GT(nbrFirstLayer[c], 0u) << name;
        EXPECT_LT(maxError[c], 1e-13 * maxValue[c]) << name;
        ++c;
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    auto const result = RUN_ALL_TESTS();
    sim_.reset();
    return result;
}
