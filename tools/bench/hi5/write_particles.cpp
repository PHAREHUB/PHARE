
#include "benchmark/benchmark.h"

#include "diagnostic/detail/h5writer.h"
#include "diagnostic/detail/h5_utils.h"
#include "diagnostic/diagnostic_manager.h"

#include "core/data/particles/particle.h"
#include "core/data/particles/particle_array.h"
#include "core/data/particles/particle_packer.h"

#include "phare/phare.h"

constexpr std::size_t dim = 1;

namespace PHARE::diagnostic
{
void do_bench(benchmark::State& state)
{
    using HiFile              = HighFive::File;
    using Packer              = core::ParticlePacker<dim>;
    using ContiguousParticles = core::ContiguousParticles<dim>;
    using ParticleArray       = core::ParticleArray<dim>;

    auto getSize = [](auto const& value) -> std::size_t {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (h5::is_array_dataset<ValueType, dim>)
            return value.size();
        else
            return 1u; /* not an array so value one of type ValueType*/
    };

    auto createDataSet_ = [&](auto& hi5, auto const& path, auto const size, auto const& value) {
        h5::createGroupsToDataSet(hi5.file_, path);
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (h5::is_array_dataset<ValueType, dim>)
            return hi5.file_.template createDataSet<typename ValueType::value_type>(
                path, HighFive::DataSpace(size));
        else
            return hi5.file_.template createDataSet<ValueType>(path, HighFive::DataSpace(size));
    };


    auto writeParticles = [](auto& datasets, auto const& particles) {
        datasets[0].write(particles.weight.data());
        datasets[1].write(particles.charge.data());
        datasets[2].write(particles.iCell.data());
        datasets[3].write(particles.delta.data());
        datasets[4].write(particles.v.data());
    };

    std::string path{"/lol/"};
    while (state.KeepRunning())
    {
        h5::HighFiveFile hi5("lol.lol", HiFile::ReadWrite | HiFile::Create | HiFile::Truncate);
        auto d = hi5.file_.createDataSet<float>("/No", 1);
        std::vector<decltype(d)> datasets;

        ContiguousParticles particles{100000};
        ParticleArray particleArray(100000);
        Packer{particleArray}.pack(particles);

        std::size_t part_idx = 0;
        core::apply(Packer::empty(), [&](auto const& arg) {
            datasets.emplace_back(createDataSet_(hi5, path + Packer::keys()[part_idx],
                                                 getSize(arg) * particles.size(), arg));
            part_idx++;
        });
        writeParticles(datasets, particles);
    }
}
BENCHMARK(do_bench)->Unit(benchmark::kMicrosecond);
} // namespace PHARE::diagnostic


int main(int argc, char* argv[])
{
    PHARE::SamraiLifeCycle samsam(argc, argv);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    return 0;
}
