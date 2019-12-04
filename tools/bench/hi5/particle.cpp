

#include "kul/log.hpp"
#include "kul/dbg.hpp"

#include "diagnostic/detail/highfive.h"
#include "diagnostic/detail/h5_utils.h"
#include "diagnostic/diagnostic_manager.h"

#include "diagnostic/samrai_lifecycle.h"

constexpr size_t dim = 1;

namespace PHARE::diagnostic
{
void bench()
{
    h5::HighFiveFile hi5("lol.lol");
    using Packer    = ParticlePacker<dim>;
    using Particles = ContiguousParticles<dim>;

    auto getSize = [&](auto const& value) -> size_t {
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (h5::is_array_dataset<ValueType, dim>)
            return value.size();
        else
            return 1u; /* not an array so value one of type ValueType*/
    };

    auto createDataSet = [&hi5](auto path, auto size, auto const& value) {
        H5Easy::detail::createGroupsToDataSet(hi5.file_, path);
        using ValueType = std::decay_t<decltype(value)>;
        if constexpr (h5::is_array_dataset<ValueType, dim>)
            return hi5.file_.createDataSet<typename ValueType::value_type>(
                path, HighFive::DataSpace(size * value.size()));
        else
            return hi5.file_.createDataSet<ValueType>(path, HighFive::DataSpace(size));
    };

    auto writeParticles = [&](auto path, auto& datasets, auto& particles) {
        KUL_DBG_FUNC_ENTER;
        datasets[0].write(particles.weight.data());
        datasets[1].write(particles.charge.data());
        datasets[2].write(particles.iCell.data());
        datasets[3].write(particles.delta.data());
        datasets[4].write(particles.v.data());
    };

    auto copyTo = [](auto& a, auto& idx, auto size, auto& v) {
        std::copy(a.begin(), a.begin() + size, v.begin() + (idx * size));
    };

    auto copyToContinguous = [&](auto& particles, auto& particleArray) {
        Packer packer(particleArray);
        size_t idx = 0;

        while (packer.hasNext())
        {
            auto next             = packer.next();
            particles.weight[idx] = std::get<0>(next);
            particles.charge[idx] = std::get<1>(next);
            copyTo(std::get<2>(next), idx, dim, particles.iCell);
            copyTo(std::get<3>(next), idx, dim, particles.delta);
            copyTo(std::get<4>(next), idx, 3, particles.v);
            idx++;
        }
    };

    auto d = hi5.file_.createDataSet<float>("/No", 1);
    std::vector<decltype(d)> datasets;

    std::string path{"/lol/"};
    Particles particles{100000};
    KLOG(INF) << particles.v[0];
    core::ParticleArray<dim> particleArray(100000);
    copyToContinguous(particles, particleArray);

    size_t part_idx = 0;
    core::apply(Packer::empty(), [&](auto const& arg) {
        datasets.emplace_back(
            createDataSet(path + Packer::keys()[part_idx], getSize(arg) * particles.size(), arg));
        part_idx++;
    });
    writeParticles(path, datasets, particles);

    KLOG(INF) << particles.v[0];
    KLOG(INF) << particles.v.back();
}
} // namespace PHARE::diagnostic

int main(int argc, char* argv[])
{
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    PHARE::diagnostic::bench();
    return 0;
}