#ifndef PHARE_CORE_UTILITIES_TIMESTAMPS_HPP
#define PHARE_CORE_UTILITIES_TIMESTAMPS_HPP

#include "core/def.hpp"

#include "initializer/data_provider.hpp"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>


namespace PHARE::core
{

enum class TimeStepType : std::uint8_t { constant, adaptive };


// Returns the dt such that start + dt, rounded to a double, equals end bit-wise.
inline double dt_between(double const start, double const end)
{
    // a step must move the clock forward
    if (!(end > start))
    {
        std::stringstream ss;
        ss << std::setprecision(17) << "dt_between: end " << end << " is not after start " << start;
        throw std::runtime_error(ss.str());
    }

    auto constexpr inf = std::numeric_limits<double>::infinity();

    // searches a dt such that start + dt == goal
    auto const solve = [&](double const goal) -> std::optional<double> {
        // first guess: when start >= goal / 2 the subtraction is exact (Sterbenz lemma)
        // and start + dt == goal at once; this covers every constant step and almost
        // every adaptive one
        auto dt = goal - start;
        // otherwise goal - start was rounded and start + dt may miss goal by an ulp;
        // each iteration corrects dt; a few suffice, 8 bounds the loop
        for (int iteration = 0; iteration < 8; ++iteration)
        {
            auto const landed = start + dt;
            if (landed == goal)
                return dt;
            // landed is close to goal (within a factor 2), so goal - landed is exact
            // (Sterbenz lemma): it is the miss, add it to dt
            auto const corrected = dt + (goal - landed);
            // if the miss is at most half an ulp of dt, adding it changes nothing:
            // move dt by one ulp towards goal instead
            dt = (corrected == dt) ? std::nextafter(dt, landed < goal ? inf : -inf) : corrected;
        }
        return std::nullopt;
    };

    // land exactly on end when some dt can
    if (auto const dt = solve(end))
        return *dt;
    // round-half-even tie: every start + dt near end falls halfway between two doubles
    // and rounds to the even one; end is odd so it is skipped, aim one ulp above it
    if (auto const dt = solve(std::nextafter(end, inf)))
        return *dt;

    std::stringstream ss;
    ss << std::setprecision(17) << "dt_between: no dt goes from " << start << " to " << end;
    throw std::runtime_error(ss.str());
}


inline bool reaches_target(double const time, double const target)
{
    return time == target
           or time == std::nextafter(target, std::numeric_limits<double>::infinity());
}


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
