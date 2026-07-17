#ifndef PHARE_CORE_UTILITIES_THREAD_POOL_HPP
#define PHARE_CORE_UTILITIES_THREAD_POOL_HPP

#include "core/logger.hpp"

#if __has_include("BS_thread_pool.hpp")
#define PHARE_HAVE_BSTP_THREAD_POOL 1
#include "BS_thread_pool.hpp"

#else
#error
#define PHARE_HAVE_BSTP_THREAD_POOL 0

#endif

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
#if PHARE_HAVE_BSTP_THREAD_POOL
        // PHARE_LOG_LINE_SS("n_pools " << n_pools);
        // PHARE_LOG_LINE_SS("threads_per_pool " << threads_per_pool);
        thread_pools.reserve(n_pools);
        for (std::uint8_t i = 0; i < n_pools; ++i)
            thread_pools.emplace_back(
                std::make_shared<::BS::thread_pool<::BS::tp::none>>(threads_per_pool));
#endif
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
#if PHARE_HAVE_BSTP_THREAD_POOL
        assert(thread_pools.size());
        thread_pools[first_ready_idx()]->detach_task(fn);
#else
        throw std::runtime_error("no impl");
#endif
    }

    void sync()
    {
#if PHARE_HAVE_BSTP_THREAD_POOL
        for (auto& p : thread_pools)
            p->wait();
#endif
    }

    auto all_finished() const
    {
#if PHARE_HAVE_BSTP_THREAD_POOL
        for (std::size_t i = 0; i < n_pools; ++i)
            if (!thread_pools[i]->wait_for(_1ms))
                return false;
#endif
        return true;
    }

    auto wait() const
    {
#if PHARE_HAVE_BSTP_THREAD_POOL
        for (auto& tp : thread_pools)
            tp->wait();
#endif
    }

    std::size_t first_ready_idx() // poll for first finished pool to reuse
    { // this is used if you want a sync point for all threads in a pool to be finished what
      // they're doing before moving on
#if PHARE_HAVE_BSTP_THREAD_POOL
        while (true)
        {
            for (; pool_idx < n_pools; ++pool_idx)
                if (thread_pools[pool_idx]->wait_for(_1ms))
                    return (++pool_idx) - 1; // :)
            pool_idx = 0;
        }
#endif
        throw std::runtime_error("no impl");
    }

    std::vector<std::mutex> mutices;
#if PHARE_HAVE_BSTP_THREAD_POOL
    std::vector<std::shared_ptr<::BS::thread_pool<::BS::tp::none>>> thread_pools;
#endif
private:
    std::size_t pool_idx = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_UTILITIES_MPI_H */
