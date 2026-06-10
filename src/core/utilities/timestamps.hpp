#ifndef PHARE_CORE_UTILITIES_TIMESTAMPS_HPP
#define PHARE_CORE_UTILITIES_TIMESTAMPS_HPP

#include "core/def.hpp"

#include "initializer/data_provider.hpp"

#include <cassert>
#include <cstdint>
#include <string>


namespace PHARE::core
{

enum class TimeStepType : std::uint8_t { constant, adaptive };


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
        assert(dt_ == new_dt);
        return last_time_ + dt_ * ++idx_;
    }

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
        if (cppdict::get_value(time_step_dict, "mode", TimeStepType::constant)
            == TimeStepType::adaptive)
            // dt_ vlaue at init is not important here: will be reset before the first step
            return std::make_unique<KahanTimeStamper>(0.);

        assert(time_step_dict.contains("value"));
        auto time_step = time_step_dict["value"].template to<double>();

        return std::make_unique<ConstantTimeStamper>(time_step);
    }
};


} // namespace PHARE::core

#endif /*PHARE_CORE_UTILITIES_TIMESTAMPS_H */
