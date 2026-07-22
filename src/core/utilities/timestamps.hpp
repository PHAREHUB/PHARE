#ifndef PHARE_CORE_UTILITIES_TIMESTAMPS_HPP
#define PHARE_CORE_UTILITIES_TIMESTAMPS_HPP

#include "core/def.hpp"

#include "initializer/data_provider.hpp"

#include <cassert>
#include <string>


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
    ConstantTimeStamper(double const& dt, double const& init_time = 0)
        : dt_{dt}
        , last_time_{init_time}
    {
    }

    double operator+=([[maybe_unused]] double const& new_dt) noexcept override
    {
        assert(dt_ == new_dt); // binary comparison - should never fail in this case
        // idx_ * dt_ stays an exact multiplication (no cumulative float drift); last_time_ is a
        // fixed seed (e.g. a restart time) added fresh each call, not accumulated iteratively.
        return last_time_ + dt_ * ++idx_;
    }

    // pure read - dt is fixed by config, nothing to (re)compute. See advance() below.
    NO_DISCARD double dt() const noexcept { return dt_; }

private:
    double dt_        = 0;
    double last_time_ = 0;
    std::size_t idx_  = 0;
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

    // pure read of the last dt cached by the last operator+=() - no computation, no mutation.
    NO_DISCARD double dt() const noexcept { return dt_; }

private:
    double dt_                 = 0;
    double last_time_          = 0;
    double error_compensation_ = 0; // Tracks the lost bits
};


struct TimeStamperFactory
{
    NO_DISCARD static std::unique_ptr<ITimeStamper> create(initializer::PHAREDict const& dict)
    {
        auto const& time_step_dict = dict["time_step"];

        // init_time stays 0: the stamper accumulates a *delta* from the start of this run;
        // Simulator::advance() adds startTime_ (== restart_time on restart) on top of it. Seeding
        // the stamper with restart_time too would add it twice.
        if (time_step_dict.contains("mode")
            && time_step_dict["mode"].template to<std::string>() == "adaptive")
            // dt_ seed is irrelevant: the first (varying) dt resets it on the first step
            return std::make_unique<KahanTimeStamper>(0.);

        assert(time_step_dict.contains("value"));
        auto time_step = time_step_dict["value"].template to<double>();

        return std::make_unique<ConstantTimeStamper>(time_step);
    }
};


} // namespace PHARE::core

#endif /*PHARE_CORE_UTILITIES_TIMESTAMPS_H */
