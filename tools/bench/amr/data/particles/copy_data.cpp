
#include "phare/phare.hpp" // samrai lifecycle

#include "amr/data/particles/particles_data.hpp"

#include "tools/bench/core/bench.hpp"

#include <cassert>

namespace PHARE::amr::bench
{
template<std::size_t dim, std::size_t interp, std::size_t op>
class ParticleDataCopy : public benchmark::Fixture
{
    auto static constexpr opts = PHARE::SimOpts<>{dim, interp};
    using PHARE_Types          = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t         = PHARE_Types::GridLayout_t;
    using ParticleArray        = PHARE_Types::ParticleArray_t;
    using ParticlesData        = PHARE::amr::ParticlesData<ParticleArray>;

public:
    void SetUp(::benchmark::State const& state) override
    {
        constexpr static std::uint32_t ppc = 1e1;

        data = std::make_unique<_DATA_>();

        assert(sourceBox == data->sourceDomain);

        data->sourceData.domainParticles = PHARE::core::bench::make_particles<dim>(ppc, sourceBox);
        assert(data->sourceData.domainParticles.size() == ppc * sourceBox.size());
        assert(data->destData.domainParticles.size() == 0);
    }

    void TearDown(::benchmark::State const& /*state*/) override {}

    void copy_data(::benchmark::State&);

private:
    static inline PHARE::amr::Box<int, dim> const sourceBox{PHARE::core::ConstArray<int, dim>(0),
                                                            PHARE::core::ConstArray<int, dim>(19)};
    static inline PHARE::amr::Box<int, dim> const destBox{PHARE::core::ConstArray<int, dim>(10),
                                                          PHARE::core::ConstArray<int, dim>(19)};


    struct _DATA_ // gbench static init doesn't play well with SAMRAI
    {
        SAMRAI::tbox::Dimension dimension{dim};
        SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(dimension)};

        SAMRAI::hier::Box sourceDomain{sourceBox};
        ParticlesData sourceData{sourceDomain, ghost, "name"};

        SAMRAI::hier::Box destDomain{destBox};
        ParticlesData destData{destDomain, ghost, "name"};
    };

    std::unique_ptr<_DATA_> data;
};

template<std::size_t dim, std::size_t interp, std::size_t op>
void ParticleDataCopy<dim, interp, op>::copy_data(::benchmark::State& state)
{
    for (auto _ : state)
    {
        this->data->destData.copy(this->data->sourceData);
        data->destData.domainParticles.clear();
    }
}

BENCHMARK_TEMPLATE_DEFINE_F(ParticleDataCopy, _2_1_1, 2, 1, 1)(benchmark::State& state)
{
    copy_data(state);
}
BENCHMARK_REGISTER_F(ParticleDataCopy, _2_1_1)->Unit(benchmark::kNanosecond);


// do the same thing twice to show variance (if any)
BENCHMARK_TEMPLATE_DEFINE_F(ParticleDataCopy, _2_1_2, 2, 1, 2)(benchmark::State& state)
{
    copy_data(state);
}
BENCHMARK_REGISTER_F(ParticleDataCopy, _2_1_2)->Unit(benchmark::kNanosecond);

} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam(argc, argv);

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
