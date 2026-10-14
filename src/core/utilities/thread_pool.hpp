#ifndef PHARE_CORE_UTILITIES_THREAD_POOL_HPP
#define PHARE_CORE_UTILITIES_THREAD_POOL_HPP

#include "BS_thread_pool.hpp"

#include <mutex>
#include <chrono>
#include <memory>
#include <vector>
#include <cassert>
#include <optional>
#include <type_traits>

namespace PHARE::core
{

struct ThreadPool
{
    static inline std::size_t n_pools          = 1;
    static inline std::size_t threads_per_pool = 1;

    auto static inline const _1ms = std::chrono::milliseconds(1);

    ThreadPool()
        : mutices(n_pools)
    {
        thread_pools.reserve(n_pools);
        for (std::size_t i = 0; i < n_pools; ++i)
            thread_pools.emplace_back(
                std::make_shared<::BS::thread_pool<::BS::tp::none>>(threads_per_pool));
    }

    static ThreadPool& INSTANCE()
    {
        static ThreadPool i;
        return i;
    }

    auto& get_pool(std::size_t const idx)
    {
        assert(n_pools > 0);
        assert(threads_per_pool > 0);
        assert(idx < thread_pools.size());
        return *thread_pools[idx];
    }

    auto& get_mutex(std::size_t const idx)
    {
        assert(idx < mutices.size());
        return mutices[idx];
    }

    static auto& pool(std::size_t const idx = 0) { return INSTANCE().get_pool(idx); }
    static auto& mutex(std::size_t const idx = 0) { return INSTANCE().get_mutex(idx); }

    void async(auto&& fn)
    {
        assert(thread_pools.size());
        thread_pools[first_ready_idx()]->detach_task(fn);
    }

    void sync()
    {
        for (auto& p : thread_pools)
            p->wait();
    }

    // patch_exec / tile_exec: view = make_view(i) is built once per patch into the slot of the
    // idle pool the patch goes to (a pool works one patch at a time, so its slot is free again
    // once it is idle). Views are only ever built and destroyed on this (calling) thread: they
    // may copy python callables (e.g. Electromag's initializer) whose copy/destruction takes the
    // GIL, which this thread holds while waiting in sync() -> deadlock if done in a task.
    // The slots are allocated once per call site, not per call.

    // fn(view) for every patch i in [0, n_patches) has run on return, a task per patch
    void patch_exec(std::size_t const n_patches, auto const& make_view, auto const& fn)
    {
        using View  = std::remove_cvref_t<decltype(make_view(std::size_t{0}))>;
        auto& slots = view_slots_<View>();

        for (std::size_t i = 0; i < n_patches; ++i)
        {
            auto const pool_idx = first_ready_idx();
            auto& view          = slots[pool_idx].emplace(make_view(i));
            get_pool(pool_idx).detach_task([&fn, &view] { fn(view); });
        }
        sync();
        for (auto& slot : slots) // don't keep views of patch data past this call
            slot.reset();
    }

    // fn(view, t) for every tile t in [0, n_tiles(view)) of every patch i has run on return, the
    // tiles of a patch are tasks on its pool's threads
    void tile_exec(std::size_t const n_patches, auto const& make_view, auto const& n_tiles,
                   auto const& fn)
    {
        using View  = std::remove_cvref_t<decltype(make_view(std::size_t{0}))>;
        auto& slots = view_slots_<View>();

        for (std::size_t i = 0; i < n_patches; ++i)
        {
            auto const pool_idx = first_ready_idx();
            auto& view          = slots[pool_idx].emplace(make_view(i));
            auto const tile_cnt = n_tiles(view);
            for (std::size_t t = 0; t < tile_cnt; ++t)
                get_pool(pool_idx).detach_task([&fn, &view, t] { fn(view, t); });
        }
        sync();
        for (auto& slot : slots) // don't keep views of patch data past this call
            slot.reset();
        sync();
    }

    auto all_finished() const
    {
        for (std::size_t i = 0; i < n_pools; ++i)
            if (!thread_pools[i]->wait_for(_1ms))
                return false;
        return true;
    }


    auto is_finished(std::size_t const idx) const
    {
        if (!thread_pools[idx]->wait_for(_1ms))
            return false;

        return true;
    }

    auto wait() const
    {
        for (auto& tp : thread_pools)
            tp->wait();
    }

    std::size_t first_ready_idx()
    {
        // poll for first finished pool to reuse
        // this is used if you want a sync point for all
        // threads in a pool to be finished what they're doing before moving on
        while (true)
        {
            for (; pool_idx < n_pools; ++pool_idx)
                if (thread_pools[pool_idx]->wait_for(_1ms))
                    return (++pool_idx) - 1; // :)
            pool_idx = 0;
        }
    }

    auto& operator()() { return thread_pools[first_ready_idx()]; }

    std::vector<std::mutex> mutices;
    std::vector<std::shared_ptr<::BS::thread_pool<::BS::tp::none>>> thread_pools;

private:
    // per View type, not per call site: calls don't nest, each empties its slots before returning
    template<typename View>
    auto& view_slots_()
    {
        static std::vector<std::optional<View>> slots(thread_pools.size());
        assert(slots.size() == thread_pools.size());
        return slots;
    }

    std::size_t pool_idx = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_UTILITIES_MPI_H */
