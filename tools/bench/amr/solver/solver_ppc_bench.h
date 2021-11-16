#ifndef PHARE_BENCH_AMR_SOLVER_PCC_BENCH_H
#define PHARE_BENCH_AMR_SOLVER_PCC_BENCH_H

#include "bench/core/bench.h"

namespace PHARE::amr::bench
{
template<typename HybridPatchState>
struct HybridPatchStateDecomposer
{
    using View  = typename HybridPatchState::view_t;
    using Batch = std::vector<View>;

    View static make_view(HybridPatchState& state, std::size_t r, std::size_t t)
    {
        auto& layout = state.layout;
        return {state};
    }


    template<typename Type>
    auto static decompose(std::vector<std::unique_ptr<HybridPatchState>>& states,
                          std::size_t threads, std::integral_constant<Type, 1> ic)
    {
        assert(false);
        return std::vector<Batch>{};
    }

    template<typename Type>
    auto static decompose(std::vector<std::unique_ptr<HybridPatchState>>& states,
                          std::size_t threads, std::integral_constant<Type, 2> ic)
    {
        std::vector<Batch> batches(threads);

        std::size_t per_thread = 0, modulo = 0;
        {
            std::size_t rows = 0;
            for (auto& state : states)
                rows += state->layout.nbrCells()[1];
            per_thread = rows / threads;
            modulo     = rows % threads;
        }

        std::size_t state_idx = 0;
        for (auto& batch : batches)
        {
            std::size_t rows_taken   = 0;
            std::int64_t rows_remain = per_thread;
            while (rows_remain > 0)
            {
                std::size_t rows_avail = states[state_idx]->layout.nbrCells()[1];
                if (rows_avail >= rows_remain)
                {
                    batch.emplace_back(make_view(*states[state_idx], rows_avail, rows_taken));
                    rows_taken += rows_avail;
                    rows_remain -= rows_avail;
                }
                else
                {
                    // todo
                }
            }
            ++state_idx;
        }

        return batches;
    }

    template<typename Type>
    auto static decompose(std::vector<std::unique_ptr<HybridPatchState>>& states,
                          std::size_t threads, std::integral_constant<Type, 3> ic)
    {
        std::vector<Batch> batches(threads);
        assert(false);
        return batches;
    }
}; // namespace PHARE::amr::bench


template<typename HybridPatchState>
auto decompose(std::vector<std::unique_ptr<HybridPatchState>>& states, std::size_t threads)
{
    return HybridPatchStateDecomposer<HybridPatchState>::decompose(
        states, threads, std::integral_constant<std::size_t, HybridPatchState::dimension>{});
}

} // namespace PHARE::amr::bench

#endif /*PHARE_BENCH_AMR_SOLVER_PCC_BENCH_H*/
