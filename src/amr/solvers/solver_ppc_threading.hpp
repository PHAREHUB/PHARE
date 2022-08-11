#ifndef PHARE_SOLVER_PPC_THREADING_HPP
#define PHARE_SOLVER_PPC_THREADING_HPP

#include <algorithm>

#include "core/utilities/types.hpp"
#include "core/utilities/thread_pool.hpp"


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


    struct PatchView
    {
        GridLayout const layout;
        typename Electromag::view_t em;
        IonsView ions;

        std::size_t const total_particles
            = core::sum_from(ions, [](auto const& pop) { return pop.domainParticles().size(); });
    };


public:
    DefaultSolverPPCThreadingStrategy(std::size_t n_threads)
        : n_threads{n_threads}
        , pool{n_threads}
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
            ion_patch_views.emplace_back(layout,                        //
                                         hybridState.electromag.view(), //
                                         IonsView::make(hybridState.ions));
        }
    }

    void solve(initializer::PHAREDict updaterDict, double dt, core::UpdaterMode mode)
    {
        auto patch_view_ptrs
            = core::generate([](auto& patch_view) { return &patch_view; }, ion_patch_views);

        if (n_threads > 1)
            std::sort(patch_view_ptrs.begin(), patch_view_ptrs.end(),
                      [](auto& a, auto& b) { return a->total_particles > b->total_particles; });

        std::vector<std::vector<PatchView*>> patch_views_per_thread(n_threads);
        {
            std::size_t thread_idx = 0;
            for (auto const& patch_view_ptr : patch_view_ptrs)
            {
                patch_views_per_thread[thread_idx++].push_back(patch_view_ptr);
                if (thread_idx == n_threads)
                    thread_idx = 0;
            }
        }

        auto thread_fn = [&](std::size_t thread_idx) {
            IonUpdater_t ionUpdater{updaterDict};
            for (auto& patch_view : patch_views_per_thread[thread_idx])
            {
                ionUpdater.updatePopulations(patch_view->ions, patch_view->em, patch_view->layout,
                                             dt, mode);
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
    std::vector<core::aggregate_adapter<PatchView>> ion_patch_views;
};

} // namespace PHARE::solver


#endif
