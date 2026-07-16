
#include "tools/bench/core/bench.hpp"
#include "core/utilities/timestamps.hpp"

#include "benchmark/benchmark.h"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>

namespace PHARE::core::bench
{
template<typename T>
std::string to_string_with_precision(T const& a_value, std::size_t const len)
{
    std::ostringstream out;
    out.precision(len);
    out << std::fixed << a_value;
    return out.str();
}

// accuracy test adapted from increment_error.cpp: repeatedly increments a
// TimeStamper by a fixed dt and reports the first step at which the
// accumulated time acquires spurious low-order digits.
template<typename Stamper>
void check_accuracy_constant_dt(std::string const& name, double const time_step = .001,
                                 std::size_t const time_step_nbr = 3000000)
{
    Stamper stamper{time_step};

    for (std::size_t i = 0; i < time_step_nbr; i++)
    {
        double const t = (stamper += time_step);
        auto const timeStamp = to_string_with_precision(t, 10);
        if (timeStamp.back() != '0')
        {
            std::cout << name << " diverged at step " << i << " : " << timeStamp << std::endl;
            return;
        }
    }
    std::cout << name << " did not diverge after " << time_step_nbr << " steps" << std::endl;
}

// A "Variable" stamper only takes its incremental (as opposed to recomputed
// from scratch) code path when new_dt actually differs from the previous dt.
// Alternating dt by +-1e-16 forces that path on every step while keeping the
// true accumulated value indistinguishable from time_step_nbr * time_step at
// the precision checked below, so any divergence still reflects the
// stamper's own floating-point handling rather than a genuinely different dt.
template<typename Stamper>
void check_accuracy_variable_dt(std::string const& name, double const time_step = .001,
                                 std::size_t const time_step_nbr = 3000000)
{
    Stamper stamper{time_step};

    for (std::size_t i = 0; i < time_step_nbr; i++)
    {
        double const dt = time_step + (i % 2 == 0 ? 1e-16 : -1e-16);
        double const t  = (stamper += dt);
        auto const timeStamp = to_string_with_precision(t, 10);
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
    // ConstantTimeStamper asserts new_dt == dt_, so it cannot take part in
    // the variable-dt check below; it is only exercised with a fixed dt.
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

// cost as used in practice, i.e. through the ITimeStamper interface returned
// by TimeStamperFactory (see simulator.hpp: (*timeStamper) += dt)
template<typename Stamper>
void step_virtual(benchmark::State& state)
{
    constexpr double time_step = .001;
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
