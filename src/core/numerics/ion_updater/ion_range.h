#ifndef PHARE_ION_RANGE_H
#define PHARE_ION_RANGE_H

#include <memory>
#include <unordered_set>

#include "core/utilities/range/ranges.h"
#include "core/utilities/range/range_replacer.h"

#include "core/logger.h"


namespace PHARE::core
{
template<typename Iterator, typename GridLayout, typename Electromag, typename PopulationView>
struct IonUpdaterRange : public Range<Iterator>
{
    using Super    = Range<Iterator>;
    using iterator = Iterator;

    IonUpdaterRange(Range<Iterator> const& range, GridLayout& layout_, Electromag& em_,
                    std::shared_ptr<PopulationView>& view_)
        : Super{range}
        , layout{&layout_}
        , em{&em_}
        , view{view_}
    {
    }

    GridLayout* layout;
    Electromag* em;
    std::shared_ptr<PopulationView> view;
};


template<typename IonUpdaterRange_t>
struct UpdaterParticleRanges
{
    UpdaterParticleRanges(IonUpdaterRange_t& domain_)
        : domain{domain_}
    {
    }

    UpdaterParticleRanges(UpdaterParticleRanges const&) = default;

    IonUpdaterRange_t domain;
    std::shared_ptr<IonUpdaterRange_t> patch_ghost, level_ghost;
};


template<typename IonUpdaterRange_t, typename ThreadRanges> //
auto pin_ghosts_to_first_range(ThreadRanges& thread_ranges)
{
    auto make_ghost_range = [](auto& ghosts, auto& ranges) -> std::shared_ptr<IonUpdaterRange_t> {
        return std::make_shared<aggregate_adapter<IonUpdaterRange_t>>(
            makeRange(ghosts), *ranges.domain.layout, *ranges.domain.em, ranges.domain.view);
    };

    std::unordered_set<void*> assigned_patch_ghost;

    for (auto& ranges_vec : thread_ranges)
        for (auto& ranges : ranges_vec)
        {
            auto& pop = *ranges.domain.view;

            if (assigned_patch_ghost.count(pop.patch_ghost) == 0)
            {
                ranges.patch_ghost = make_ghost_range(*pop.patch_ghost, ranges);
                ranges.level_ghost = make_ghost_range(*pop.level_ghost, ranges);
                assigned_patch_ghost.emplace(pop.patch_ghost);
            }
        }
}

} // namespace PHARE::core

namespace PHARE::core::detail
{
template<typename IonUpdaterRange_t, typename GridLayout, typename Electromag, typename IonPopView>
auto make_ranges(
    std::vector<std::tuple<GridLayout, Electromag, std::vector<std::shared_ptr<IonPopView>>>>&
        views_per_patch,
    std::size_t n_threads, std::size_t n_pops)
{
    abort_if(views_per_patch.size() == 0 || n_pops == 0);

    std::vector<std::vector<std::tuple<GridLayout*, Electromag*, std::shared_ptr<IonPopView>>>>
        views_per_pop(n_pops);

    for (auto& [layout, em, pops] : views_per_patch)
        for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx)
            views_per_pop[pop_idx].emplace_back(&layout, &em, pops[pop_idx]);

    return generate(
        [&](auto i) {
            auto particles = generate(
                [&](auto& tuple) {
                    auto& [layout, em, view] = tuple;
                    return std::forward_as_tuple(layout, em, view, view->domain);
                },
                views_per_pop[i]);
            auto accesor = [&](auto& tuple) -> auto& { return *std::get<3>(tuple); };
            auto builder = [&](auto&& range, auto& tuple) {
                auto& [layout, em, view, _] = tuple;
                return IonUpdaterRange_t{range, *layout, *em, view};
            };
            return PHARE::core::make_balanced_ranges(particles, n_threads, accesor, builder);
        },
        views_per_pop.size());
}

}

namespace PHARE::core
{
template<typename GridLayout, typename Electromag, typename IonPopView> // patches // populations
auto updater_ranges_per_thread(
    std::vector<std::tuple<GridLayout, Electromag, std::vector<std::shared_ptr<IonPopView>>>>&
        views,
    std::size_t n_threads = 1)
{
    // resetMoments(ions);  could go here

    using ParticleArray = typename IonPopView::particle_array_type;
    using IonUpdaterRange_t
        = IonUpdaterRange<typename ParticleArray::iterator, GridLayout, Electromag, IonPopView>;
    using UpdaterRange_t = UpdaterParticleRanges<IonUpdaterRange_t>;
    using ThreadViews    = std::vector<std::vector<UpdaterRange_t>>;

    ThreadViews thread_views(n_threads);
    if (views.size() == 0)
        return thread_views;

    // assumes all patches have the same number of populations
    auto n_pops                = std::get<2>(views[0]).size();
    auto thread_ranges_per_pop = detail::make_ranges<IonUpdaterRange_t>(views, n_threads, n_pops);

    for (std::size_t t_idx = 0; t_idx < n_threads; ++t_idx)
        for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx)
            for (auto& range : thread_ranges_per_pop[pop_idx][t_idx])
                thread_views[t_idx].emplace_back(range);

    pin_ghosts_to_first_range<IonUpdaterRange_t>(thread_views);

    return thread_views;
}


} // namespace PHARE::core


#endif // ION_UPDATER_H
