#ifndef PHARE_AMR_LOAD_BALANCER_LOAD_BALANCER_DETAILS_HPP
#define PHARE_AMR_LOAD_BALANCER_LOAD_BALANCER_DETAILS_HPP

#include <string>
#include <cstdint>

#include "core/logger.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::amr
{
struct LoadBalancerDetailsDefaults
{
    // is instantiated in constexpr context so following fields are constexpr in that context
    std::size_t next_rebalance_backoff_multiplier = 2;
    std::size_t next_rebalance                    = 200;
    std::size_t max_next_rebalance                = 1000;
    double tolerance                              = .05;
};

struct LoadBalancerDetails
{
    LoadBalancerDetailsDefaults constexpr static defaults{};

    bool const active    = false;
    bool const automatic = false;
    bool const on_init   = false;

    std::size_t const every = 0;
    std::string const mode;

    double const tolerance = defaults.tolerance;

    std::size_t next_rebalance_backoff_multiplier = defaults.next_rebalance_backoff_multiplier;
    std::size_t next_rebalance                    = defaults.next_rebalance;
    std::size_t max_next_rebalance                = defaults.max_next_rebalance;

    LoadBalancerDetails static FROM(initializer::PHAREDict const& dict)
    {
        return {
            cppdict::get_value(dict, "active", false),
            cppdict::get_value(dict, "auto", false),
            cppdict::get_value(dict, "on_init", false),
            cppdict::get_value(dict, "every", std::size_t{0}),
            cppdict::get_value(dict, "mode", std::string{"nppc"}),
            cppdict::get_value(dict, "tolerance", defaults.tolerance),
            cppdict::get_value(dict, "next_rebalance_backoff_multiplier",
                               defaults.next_rebalance_backoff_multiplier),
            cppdict::get_value(dict, "next_rebalance", defaults.next_rebalance),
            cppdict::get_value(dict, "max_next_rebalance", defaults.max_next_rebalance),
        };
    }
};

} // namespace PHARE::amr

#endif /* PHARE_AMR_LOAD_BALANCER_LOAD_BALANCER_DETAILS_HPP */
