#ifndef PHARE_CORE_UTILITIES_THREAD_POOL_HPP
#define PHARE_CORE_UTILITIES_THREAD_POOL_HPP

#include "BS_thread_pool.hpp"

#include <mutex>
#include <chrono>
#include <memory>
#include <vector>
#include <cassert>
#include <utility>
#include <optional>
#include <exception>
#include <stdexcept>
#include <type_traits>

namespace PHARE::core
{

struct ThreadPool
{
    static inline std::size_t n_pools          = 1;
    static inline std::size_t threads_per_pool = 1;

    auto static inline const _1ms = std::chrono::milliseconds(1);

    ThreadPool() { configure(n_pools, threads_per_pool); }

    // (re)builds the pools, only while no task is in flight (e.g. Simulator::setup)
    void configure(std::size_t const pools, std::size_t const threads)
    {
        if (pools == 0 or threads == 0)
            throw std::runtime_error("ThreadPool: pools and threads_per_pool must be > 0");

        stop();
        n_pools          = pools;
        threads_per_pool = threads;
        mutices          = std::vector<std::mutex>(n_pools);
        thread_pools.reserve(n_pools);
        for (std::size_t i = 0; i < n_pools; ++i)
            thread_pools.emplace_back(
                std::make_shared<::BS::thread_pool<::BS::tp::none>>(threads_per_pool));
    }

    // joins and destroys every pool's threads, configure() is needed before further use
    void stop()
    {
        thread_pools.clear(); // BS's destructor waits for queued tasks, then joins
        pool_idx = 0;
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
        thread_pools[first_ready_idx()]->detach_task(guard(fn));
    }

    // BS's worker loop swallows exceptions: every task goes through guard(fn), which keeps the
    // first one thrown for sync() to rethrow on the calling thread, where it takes the same path
    // as in a serial run (PHARE_LOG_ERROR -> MPI_Abort or emergency dump)
    auto guard(auto&& fn)
    {
        return [this, fn = std::forward<decltype(fn)>(fn)]() mutable {
            try
            {
                fn();
            }
            catch (...)
            {
                std::scoped_lock lock{error_mutex_};
                if (!error_)
                    error_ = std::current_exception();
            }
        };
    }

    // all tasks of all pools have run on return, rethrows the first exception of any of them
    void sync()
    {
        join_();
        rethrow_();
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
            get_pool(pool_idx).detach_task(guard([&fn, &view] { fn(view); }));
        }
        join_();
        for (auto& slot : slots) // don't keep views of patch data past this call
            slot.reset();
        rethrow_();
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
                get_pool(pool_idx).detach_task(guard([&fn, &view, t] { fn(view, t); }));
        }
        join_();
        for (auto& slot : slots) // don't keep views of patch data past this call
            slot.reset();
        rethrow_();
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

    void wait() { sync(); }

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
        static std::vector<std::optional<View>> slots;
        if (slots.size() != thread_pools.size()) // pools were reconfigured, slots are empty
            slots.resize(thread_pools.size());
        return slots;
    }

    void join_()
    {
        for (auto& p : thread_pools)
            p->wait();
    }

    void rethrow_()
    {
        if (auto error = std::exchange(error_, nullptr))
            std::rethrow_exception(error);
    }

    std::size_t pool_idx = 0;
    std::mutex error_mutex_;
    std::exception_ptr error_ = nullptr;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_UTILITIES_MPI_H */
