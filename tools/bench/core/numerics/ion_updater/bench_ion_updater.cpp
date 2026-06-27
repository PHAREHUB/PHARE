
#ifndef PHARE_UPDATER_IMPL
#error // no
#endif


#include "bench_updater.hpp"

namespace PHARE::core
{

BENCHMARK_TEMPLATE_DEFINE_F(UpdaterBencherCopy, copy_3_1_x,
                            Params<3, 1, PHARE_UPDATER_IMPL>)(benchmark::State& st)
{
    (*this)(st);
}
BENCHMARK_REGISTER_F(UpdaterBencherCopy, copy_3_1_x)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE_DEFINE_F(UpdaterBencherRef, ref_3_1_x,
                            Params<3, 1, PHARE_UPDATER_IMPL>)(benchmark::State& st)
{
    (*this)(st);
}
BENCHMARK_REGISTER_F(UpdaterBencherRef, ref_3_1_x)->Unit(benchmark::kMicrosecond);

} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
