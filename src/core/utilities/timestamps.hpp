#ifndef PHARE_CORE_UTILITIES_TIMESTAMPS_HPP
#define PHARE_CORE_UTILITIES_TIMESTAMPS_HPP

#include <string>
#include <cassert>
#include <cstdint>

#include "core/logger.hpp"
#include "initializer/data_provider.hpp"
#include "core/def.hpp"

namespace PHARE::core
{
struct ITimeStamper
{
    virtual double operator+=(double const& new_dt) noexcept = 0;

    virtual ~ITimeStamper() {}
};


class ConstantTimeStamper : public ITimeStamper
{
public:
    ConstantTimeStamper(double const& dt, std::size_t const& init_idx = 0)
        : dt_{dt}
        , idx_{init_idx}
    {
    }

    double operator+=([[maybe_unused]] double const& new_dt) noexcept override
    {
        assert(dt_ == new_dt); // binary comparison - should never fail in this case
        return dt_ * ++idx_;
    }

private:
    double dt_       = 0;
    std::size_t idx_ = 0;
};

// Accumulates the (possibly varying) dt actually used at each step, for adaptive time stepping.
class VariableTimeStamper : public ITimeStamper
{
public:
    VariableTimeStamper(double const& dt, double const& init_time = 0)
        : dt_{dt}
        , last_change_{init_time}
        , last_time_(init_time)
    {
    }

    double operator+=(double const& new_dt) noexcept override
    {
        assert(new_dt > 0);

        if (new_dt != dt_) // not sure if safe, possibly
        {
            dt_          = new_dt;
            last_change_ = last_time_;
            n_same_      = 0;
        }
        return (last_time_ = last_change_ + (dt_ * ++n_same_));
    }

private:
    std::size_t n_same_ = 0;
    double dt_          = 0;
    double last_change_ = 0;
    double last_time_   = 0;
};

struct TimeStamperFactory
{
    NO_DISCARD static std::unique_ptr<ITimeStamper> create(initializer::PHAREDict const& dict)
    {
        auto const& time_step_dict = dict["time_step"];
        if (time_step_dict.contains("mode")
            && time_step_dict["mode"].template to<std::string>() == "adaptive")
            // dt_ seed is irrelevant: the first (varying) dt resets it on the first step
            return std::make_unique<VariableTimeStamper>(0.);

        assert(time_step_dict.contains("value"));
        auto time_step  = time_step_dict["value"].template to<double>();
        std::size_t idx = 0;

        return std::make_unique<ConstantTimeStamper>(time_step, idx);
    }
};


} // namespace PHARE::core

#endif /*PHARE_CORE_UTILITIES_TIMESTAMPS_H */
