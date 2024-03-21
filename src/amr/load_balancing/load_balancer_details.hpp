#ifndef PHARE_AMR_LOAD_BALANCER_LOAD_BALANCER_DETAILS_HPP
#define PHARE_AMR_LOAD_BALANCER_LOAD_BALANCER_DETAILS_HPP

#include <string>
#include <cstdint>

#include "core/logger.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::amr
{
struct LoadBalancerDetails
{
    bool const active    = false;
    bool const automatic = false;
    bool const on_init   = false;

    std::size_t const every = 0;
    std::string const mode;

    double const tolerance = .05;

    LoadBalancerDetails static FROM(initializer::PHAREDict const& dict)
    {
        return {cppdict::get_value(dict, "active", false),
                cppdict::get_value(dict, "auto", false),
                cppdict::get_value(dict, "on_init", false),
                cppdict::get_value(dict, "every", std::size_t{0}),
                cppdict::get_value(dict, "mode", std::string{"nppc"}),
                cppdict::get_value(dict, "tol", 0.05)};
    }
};

} // namespace PHARE::amr

#endif /* PHARE_AMR_LOAD_BALANCER_LOAD_BALANCER_DETAILS_HPP */
