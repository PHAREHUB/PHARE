#ifndef PHARE_SOLVER_PPC_THREADING_HPP
#define PHARE_SOLVER_PPC_THREADING_HPP

#include <thread>

#include "core/utilities/thread_pool.hpp"
#include "core/utilities/range/ranges.hpp"


namespace PHARE::solver
{
template<typename HybridModel>
class DefaultSolverPPCThreadingStrategy
{
    using This = DefaultSolverPPCThreadingStrategy;

    using Electromag       = typename HybridModel::electromag_type;
    using Ions             = typename HybridModel::ions_type;
    using ParticleArray    = typename Ions::particle_array_type;
    using VecFieldT        = typename HybridModel::vecfield_type;
    using GridLayout       = typename HybridModel::gridlayout_type;
    using ResourcesManager = typename HybridModel::resources_manager_type;

    using IonsView     = core::IonsView<ParticleArray, VecFieldT, GridLayout>;
    using IonPopView   = core::IonPopulationView<ParticleArray, VecFieldT, GridLayout>;
    using IonUpdater_t = core::IonUpdater<IonsView, typename Electromag::view_t, GridLayout>;

    using PatchView = std::tuple<GridLayout, typename Electromag::view_t, IonsView>;


    struct IonView
    {
        IonPopView& pop;

        auto from(std::vector<std::shared_ptr<IonPopView>>& ion_pop_views)
        {
            return core::generate_from([](auto& pop_view) { return IonView(*pop_view); },
                                       ion_pop_views);
        }
    };

public:
    DefaultSolverPPCThreadingStrategy(std::size_t n_threads)
        : n_threads{n_threads}
        , pool{n_threads - 1}
    {
    }


    template<typename Level>
    void build_from(Level& level, HybridModel& hybridModel, ResourcesManager& resourcesManager)
    {
        ion_patch_views.clear();
        auto& hybridState = hybridModel.state;
        for (auto& patch : *level)
        {
            auto _      = resourcesManager.setOnPatch(*patch, hybridState);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            ion_patch_views.emplace_back(layout, hybridState.electromag.view(),
                                         IonsView::make(hybridState.ions));
        }
    }

    void solve(initializer::PHAREDict updaterDict, double dt, core::UpdaterMode mode)
    {
        auto patch_views
            = core::generate_from([](auto& patch_view) { return &patch_view; }, ion_patch_views);

        std::vector<std::tuple<PatchView*, ParticleArray*>> domain_particles;
        for (auto const& patch_view : patch_views)
            for (auto const& pop : std::get<2>(*patch_view).pops)
                domain_particles.push_back(std::make_tuple(patch_view, &pop.domainParticles()));

        auto ranges = make_balanced_ranges(
            domain_particles, n_threads, [](auto& el) -> auto& { return *std::get<1>(el); },
            [](auto&& range, auto& el) { return std::make_tuple(range, el); });

        auto thread_fn = [&](std::size_t thread_idx) {
            IonUpdater_t ionUpdater{updaterDict};
            for (auto& range_tuple : ranges[thread_idx])
            {
                auto& [range, patch_view_tuple] = range_tuple;
                auto& [patch_view, _]           = patch_view_tuple;
                auto& [layout, electromag, __]  = *patch_view;
                // ionUpdater.updatePopulations(range, electromag, layout, dt, mode); // fix
            }
        };

        for (std::size_t i = 1; i < n_threads; ++i)
            pool.submit([&, i]() { thread_fn(i); });

        thread_fn(0);
        pool.wait_for_tasks();
    }

private:
    std::size_t n_threads = 1;
    thread_pool pool;
    std::vector<PatchView> ion_patch_views;
};

} // namespace PHARE::solver


#endif
