#ifndef PHARE_CORE_UTILITIES_TIMESTAMPS_HPP
#define PHARE_CORE_UTILITIES_TIMESTAMPS_HPP

#include "core/def.hpp"
#include "core/logger.hpp"
#include "initializer/data_provider.hpp"

#include <string>
#include <cassert>
#include <cstdint>

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

// No error compensation: reduces to plain incremental addition whenever dt
// actually changes step to step, and is therefore subject to the same
// accumulated rounding drift as naive floating-point summation.
class NaiveTimeStamper : public ITimeStamper
{
public:
    NaiveTimeStamper(double const& dt, double const& init_time = 0)
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

class KahanTimeStamper : public ITimeStamper
{
public:
    KahanTimeStamper(double const& dt, double const& init_time = 0)
        : dt_{dt}
        , last_time_(init_time)
        , error_compensation_(0.0)
    {
    }

    double operator+=(double const& new_dt) noexcept override
    {
        assert(new_dt > 0);
        dt_ = new_dt;

        // Kahan Summation Algorithm
        // 1. Subtract the accumulated floating-point error from the new increment
        double y = dt_ - error_compensation_;

        // 2. Add the corrected increment to the running total.
        // High-order bits are updated; low-order bits of 'y' may be lost here.
        double temp = last_time_ + y;

        // 3. Recover the exact low-order bits that were dropped in the addition
        error_compensation_ = (temp - last_time_) - y;

        // 4. Update the state
        return (last_time_ = temp);
    }

private:
    double dt_                 = 0;
    double last_time_          = 0;
    double error_compensation_ = 0; // Tracks the lost bits
};

// Same compensated-summation idea as KahanTimeStamper, restructured around a
// tail carried forward into the next increment rather than an error term
// subtracted from it.
class CascadedTimeStamper : public ITimeStamper
{
public:
    CascadedTimeStamper(double const& dt, double const& init_time = 0)
        : dt_{dt}
        , time_head_(init_time)
        , time_tail_(0.0)
    {
    }

    double operator+=(double const& new_dt) noexcept override
    {
        assert(new_dt > 0);
        dt_ = new_dt;

        // Cascade the new timestep into the lower-precision tail first
        double interim_tail = time_tail_ + dt_;

        // Advance the head by the combined tail values
        double next_head = time_head_ + interim_tail;

        // Calculate what failed to transfer to the head, preserving it in the tail
        time_tail_ = interim_tail - (next_head - time_head_);
        time_head_ = next_head;

        return time_head_ + time_tail_;
    }

private:
    double dt_        = 0;
    double time_head_ = 0; // Stores the macroscopic time
    double time_tail_ = 0; // Stores the sub-CFL micro-fractions
};


struct TimeStamperFactory
{
    std::string constexpr static default_time_stamper = "constant";

    NO_DISCARD static std::unique_ptr<ITimeStamper> create(initializer::PHAREDict const& dict)
    {
        assert(dict.contains("time_step"));
        auto time_step = dict["time_step"].template to<double>();

        auto type = cppdict::get_value(dict, "time_stamper", default_time_stamper);

        if (type == "constant")
            return std::make_unique<ConstantTimeStamper>(time_step);

        if (type == "naive")
            return std::make_unique<NaiveTimeStamper>(time_step);
        if (type == "kahan")
            return std::make_unique<KahanTimeStamper>(time_step);
        if (type == "cascaded")
            return std::make_unique<CascadedTimeStamper>(time_step);

        throw std::runtime_error("No TimeStamper exists for key: " + type);
    }
};


} // namespace PHARE::core

#endif /*PHARE_CORE_UTILITIES_TIMESTAMPS_H */
