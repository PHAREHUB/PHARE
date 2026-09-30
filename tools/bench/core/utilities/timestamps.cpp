
#include "core/utilities/timestamps.hpp"
#include "core/utilities/types.hpp"

#include "benchmark/benchmark.h"

#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

namespace PHARE::core::bench
{

// No error compensation: reduces to plain incremental addition whenever dt
// actually changes step to step, and is therefore subject to the same
// accumulated rounding drift as naive floating-point summation.
class NaiveTimeStamper : public ITimeStamper
{
public:
    NaiveTimeStamper(double const& dt, double const& init_time = 0)
        : dt_{dt}
        , last_change_{init_time}
        , last_time_{init_time}
    {
    }

    double operator+=(double const& new_dt) noexcept override
    {
        assert(new_dt > 0);

        if (new_dt != dt_)
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

// Same compensated-summation idea as KahanTimeStamper, restructured around a
// tail carried forward into the next increment rather than an error term
// subtracted from it.
class CascadedTimeStamper : public ITimeStamper
{
public:
    CascadedTimeStamper(double const& dt, double const& init_time = 0)
        : dt_{dt}
        , time_head_{init_time}
        , time_tail_{0.}
    {
    }

    double operator+=(double const& new_dt) noexcept override
    {
        assert(new_dt > 0);
        dt_ = new_dt;

        // cascade the new timestep into the lower-precision tail first
        double const interim_tail = time_tail_ + dt_;

        // advance the head by the combined tail values
        double const next_head = time_head_ + interim_tail;

        // what failed to transfer to the head is preserved in the tail
        time_tail_ = interim_tail - (next_head - time_head_);
        time_head_ = next_head;

        return time_head_ + time_tail_;
    }

private:
    double dt_        = 0;
    double time_head_ = 0; // macroscopic time
    double time_tail_ = 0; // sub-CFL micro-fractions
};


// accuracy test adapted from increment_error.cpp: repeatedly increments a
// TimeStamper by a fixed dt and reports the first step at which the accumulated
// time acquires spurious low-order digits.
// https://github.com/PHARCHIVE/test_snippets/blob/main/numeric/double/increment_error.cpp
template<typename Stamper>
void check_accuracy_constant_dt(std::string const& name, double const time_step = .001,
                                std::size_t const time_step_nbr = 3000000)
{
    Stamper stamper{time_step};

    for (std::size_t i = 0; i < time_step_nbr; ++i)
    {
        auto const timeStamp = to_string_with_precision(stamper += time_step, 10);
        if (timeStamp.back() != '0')
        {
            std::cout << name << " diverged at step " << i << " : " << timeStamp << std::endl;
            return;
        }
    }
    std::cout << name << " did not diverge after " << time_step_nbr << " steps" << std::endl;
}

// A stamper only takes its incremental (as opposed to recomputed from scratch)
// code path when new_dt actually differs from the previous dt. Alternating dt by
// +-1e-16 forces that path on every step while keeping the true accumulated
// value indistinguishable from time_step_nbr * time_step at the precision
// checked below, so any divergence still reflects the stamper's own
// floating-point handling rather than a genuinely different dt.
template<typename Stamper>
void check_accuracy_variable_dt(std::string const& name, double const time_step = .001,
                                std::size_t const time_step_nbr = 3000000)
{
    Stamper stamper{time_step};

    for (std::size_t i = 0; i < time_step_nbr; ++i)
    {
        double const dt      = time_step + (i % 2 == 0 ? 1e-16 : -1e-16);
        auto const timeStamp = to_string_with_precision(stamper += dt, 10);
        if (timeStamp.back() != '0')
        {
            std::cout << name << " (variable dt) diverged at step " << i << " : " << timeStamp
                      << std::endl;
            return;
        }
    }
    std::cout << name << " (variable dt) did not diverge after " << time_step_nbr << " steps"
              << std::endl;
}

void check_all_accuracy()
{
    // ConstantTimeStamper asserts new_dt == dt_, so it cannot take part in the
    // variable-dt check below; it is only exercised with a fixed dt.
    check_accuracy_constant_dt<ConstantTimeStamper>("ConstantTimeStamper");
    check_accuracy_constant_dt<NaiveTimeStamper>("NaiveTimeStamper");
    check_accuracy_constant_dt<KahanTimeStamper>("KahanTimeStamper");
    check_accuracy_constant_dt<CascadedTimeStamper>("CascadedTimeStamper");

    check_accuracy_variable_dt<NaiveTimeStamper>("NaiveTimeStamper");
    check_accuracy_variable_dt<KahanTimeStamper>("KahanTimeStamper");
    check_accuracy_variable_dt<CascadedTimeStamper>("CascadedTimeStamper");
}


// cost of incrementing the concrete stamper type directly (no virtual dispatch)
template<typename Stamper>
void step(benchmark::State& state)
{
    constexpr double time_step = .001;
    Stamper stamper{time_step};

    while (state.KeepRunning())
        benchmark::DoNotOptimize(stamper += time_step);
}

// cost as used in practice, i.e. through the ITimeStamper interface returned by
// TimeStamperFactory (see simulator.hpp: (*timeStamper_) += dt)
template<typename Stamper>
void step_virtual(benchmark::State& state)
{
    constexpr double time_step            = .001;
    std::unique_ptr<ITimeStamper> stamper = std::make_unique<Stamper>(time_step);

    while (state.KeepRunning())
        benchmark::DoNotOptimize((*stamper) += time_step);
}

BENCHMARK_TEMPLATE(step, ConstantTimeStamper)->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(step, NaiveTimeStamper)->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(step, KahanTimeStamper)->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(step, CascadedTimeStamper)->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(step_virtual, ConstantTimeStamper)->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(step_virtual, NaiveTimeStamper)->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(step_virtual, KahanTimeStamper)->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(step_virtual, CascadedTimeStamper)->Unit(benchmark::kNanosecond);

} // namespace PHARE::core::bench


int main(int argc, char** argv)
{
    PHARE::core::bench::check_all_accuracy();

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
