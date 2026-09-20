#ifndef PHARE_CORE_UTILITIES_THREAD_POOL_HPP
#define PHARE_CORE_UTILITIES_THREAD_POOL_HPP

#include "BS_thread_pool.hpp"

#include <mutex>
#include <chrono>
#include <memory>
#include <cassert>

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
    std::size_t pool_idx = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_UTILITIES_MPI_H */
