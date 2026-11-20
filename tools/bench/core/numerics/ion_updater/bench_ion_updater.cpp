#include "bench_updater.hpp"

namespace PHARE::core
{

BENCHMARK_TEMPLATE_DEFINE_F(UpdaterBencherCopy, copy_3_1, Params<3, 1>)(benchmark::State& st)
{
    (*this)(st);
}
BENCHMARK_REGISTER_F(UpdaterBencherCopy, copy_3_1)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE_DEFINE_F(UpdaterBencherCopy, copy_3_1_reference,
                            Params<3, 1, /*reference*/ true>)(benchmark::State& st)
{
    (*this)(st);
}
BENCHMARK_REGISTER_F(UpdaterBencherCopy, copy_3_1_reference)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE_DEFINE_F(UpdaterBencherInPlace, in_place_3_1, Params<3, 1>)(benchmark::State& st)
{
    (*this)(st);
}
BENCHMARK_REGISTER_F(UpdaterBencherInPlace, in_place_3_1)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE_DEFINE_F(UpdaterBencherInPlace, in_place_3_1_reference,
                            Params<3, 1, /*reference*/ true>)(benchmark::State& st)
{
    (*this)(st);
}
BENCHMARK_REGISTER_F(UpdaterBencherInPlace, in_place_3_1_reference)->Unit(benchmark::kMicrosecond);

} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
