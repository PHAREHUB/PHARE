
#include "tools/bench/core/bench.hpp"

#include "amr/data/particles/particles_data.hpp"

#include "amr/utilities/box/amr_box.hpp"

#include "phare/phare.hpp" // samrai lifecycle

namespace PHARE::amr::bench
{
template<std::size_t dim, std::size_t interp, std::size_t op>
class ParticleDataCopy : public benchmark::Fixture
{
    using PHARE_Types   = PHARE::core::PHARE_Types<dim, interp>;
    using GridLayout_t  = typename PHARE_Types::GridLayout_t;
    using ParticleArray = typename PHARE_Types::ParticleArray_t;
    using ParticlesData = PHARE::amr::ParticlesData<ParticleArray>;

public:
    void SetUp(::benchmark::State const& state) override
    {
        constexpr static std::uint32_t ppc = 1e1;

        data = std::make_unique<_DATA_>();

        assert(sourceBox == data->sourceDomain);
        assert(sourceBox == (PHARE::amr::Box<int, dim>{data->sourceDomain}));

        data->sourceData.domainParticles = PHARE::core::bench::make_particles<dim>(ppc, sourceBox);
        assert(data->sourceData.domainParticles.size() == ppc * sourceBox.size());
    }

    void TearDown(::benchmark::State const& /*state*/) override {}

    void copy_data(::benchmark::State&);

private:
    static const inline PHARE::amr::Box<int, dim> sourceBox{PHARE::core::ConstArray<int, dim>(0),
                                                            PHARE::core::ConstArray<int, dim>(999)};
    static const inline PHARE::amr::Box<int, dim> destBox{PHARE::core::ConstArray<int, dim>(500),
                                                          PHARE::core::ConstArray<int, dim>(1499)};


    struct _DATA_ // gbench static init doesn't play well with SAMRAI
    {
        SAMRAI::tbox::Dimension dimension{dim};
        SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(dimension)};

        SAMRAI::hier::Box sourceDomain{sourceBox};
        ParticlesData sourceData{sourceDomain, ghost};

        SAMRAI::hier::Box destDomain{destBox};
        ParticlesData destData{destDomain, ghost};
    };

    std::unique_ptr<_DATA_> data;
};

template<std::size_t dim, std::size_t interp, std::size_t op>
void ParticleDataCopy<dim, interp, op>::copy_data(::benchmark::State& state)
{
    for (auto _ : state)
        this->data->destData.copy(this->data->sourceData);
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
