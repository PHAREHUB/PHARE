#ifndef PHARE_AMR_LOAD_BALANCER_COARSEST_REBALANCE_CHECK_HPP
#define PHARE_AMR_LOAD_BALANCER_COARSEST_REBALANCE_CHECK_HPP

#include "amr/load_balancing/load_balancer_details.hpp"

#include <vector>
#include <cstddef>
#include <numeric>
#include <algorithm>

namespace PHARE::amr
{
/**
 * Decides, once per time step, whether L0 is explicitly rebalanced in automatic mode.
 *
 * The workloads are only gathered when a check is due, every rank holds the same state so
 * they all agree on when that is, which matters as gathering them is collective.
 *
 * After a check that finds L0 balanced the interval to the next check grows by
 * next_rebalance_backoff_multiplier up to max_next_rebalance, after a rebalance it
 * restarts from next_rebalance.
 */
class CoarsestRebalanceCheck
{
public:
    CoarsestRebalanceCheck(LoadBalancerDetails const& info)
        : info_{info}
        , interval_{at_least_one(info.next_rebalance)}
        , steps_to_next_check_{info.on_init ? 0 : interval_}
    {
    }

    // same definition as SAMRAI's flexible_load_tolerance: the most loaded rank may carry up
    // to tolerance above the average load
    bool static is_imbalanced(std::vector<double> const& workloads, double const tolerance)
    {
        if (workloads.empty())
            return false;

        auto const total = std::accumulate(workloads.begin(), workloads.end(), 0.0);
        if (total <= 0)
            return false;

        auto const average = total / workloads.size();
        auto const max     = *std::max_element(workloads.begin(), workloads.end());
        return max / average - 1 > tolerance;
    }

    template<typename GatherWorkloads>
    bool operator()(GatherWorkloads&& gather_workloads)
    {
        if (steps_to_next_check_ > 0)
        {
            --steps_to_next_check_;
            return false;
        }

        bool const rebalance = is_imbalanced(gather_workloads(), info_.tolerance);

        if (rebalance)
            interval_ = at_least_one(info_.next_rebalance);

        steps_to_next_check_ = interval_ - 1;

        if (!rebalance)
            interval_ = std::min(interval_ * at_least_one(info_.next_rebalance_backoff_multiplier),
                                 at_least_one(info_.max_next_rebalance));

        return rebalance;
    }

private:
    std::size_t static at_least_one(std::size_t const steps)
    {
        return std::max(steps, std::size_t{1});
    }

    LoadBalancerDetails const info_;
    std::size_t interval_            = 1;
    std::size_t steps_to_next_check_ = 0;
};

} // namespace PHARE::amr

#endif /* PHARE_AMR_LOAD_BALANCER_COARSEST_REBALANCE_CHECK_HPP */
