#ifndef PHARE_TEST_CORE_DATA_ION_POPULATIONS_ION_POPULATION_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_ION_POPULATIONS_ION_POPULATION_FIXTURES_HPP


#include "phare_core.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"


#include <cassert>


namespace PHARE::core
{

template<typename ParticleArray_, std::size_t interp_>
struct UsableIonsDefaultTypes
{
public:
    auto static constexpr dim    = ParticleArray_::dimension;
    auto static constexpr interp = interp_;
    SimOpts static constexpr opts{dim, interp_};

    using PHARE_Types         = PHARE::core::PHARE_Types<opts>;
    using ParticleArray_t     = ParticleArray_;
    using Array_t             = NdArrayVector<dim, double, /*c_ordering*/ true>;
    using Grid_t              = Grid<Array_t, HybridQuantity::Scalar>;
    using UsableVecField_t    = UsableVecField<dim>;
    using UsableTensorField_t = UsableTensorField<dim, /*rank = */ 2>;
    using GridLayout_t        = PHARE_Types::GridLayout_t;
    using VecField_t          = UsableVecField_t::Super;
    using TensorField_t       = UsableTensorField_t::Super;

    using IonPopulation_t = IonPopulation<ParticleArray_t, VecField_t, TensorField_t>;
};



auto inline pop_dict(std::string const& name, std::size_t const ppc = 0)
{
    initializer::PHAREDict popdict;
    popdict["name"]                           = name;
    popdict["mass"]                           = 1.0;
    auto particle_initializer                 = initializer::PHAREDict{};
    particle_initializer["nbr_part_per_cell"] = static_cast<int>(ppc);
    particle_initializer["charge"]            = 1.;
    particle_initializer["basis"]             = std::string{"cartesian"};
    popdict["particle_initializer"]           = particle_initializer;
    return popdict;
}


template<typename _defaults>
class UsableIonsPopulation : public _defaults::IonPopulation_t
{
    using GridLayout_t    = _defaults::GridLayout_t;
    using VecField_t      = _defaults::VecField_t;
    using TensorField_t   = _defaults::TensorField_t;
    using ParticleArray_t = _defaults::ParticleArray_t;
    using Super           = IonPopulation<ParticleArray_t, VecField_t, TensorField_t>;

    void set()
    {
        auto&& [_F, _M, _pd, _cd, _particles] = Super::getCompileTimeResourcesViewList();
        F.set_on(_F);
        M.set_on(_M);
        _pd.setBuffer(&particleDensity);
        _cd.setBuffer(&chargeDensity);
        _particles.setBuffer(&particles.pack());
    }

public:
    UsableIonsPopulation(initializer::PHAREDict const& dict, GridLayout_t const& layout)
        : Super{dict}
        , particleDensity{this->name() + "_particleDensity", layout, HybridQuantity::Scalar::rho}
        , chargeDensity{this->name() + "_chargeDensity", layout, HybridQuantity::Scalar::rho}
        , F{this->name() + "_flux", layout, HybridQuantity::Vector::V}
        , M{this->name() + "_momentumTensor", layout, HybridQuantity::Tensor::M}
        , particles{this->name(), layout.AMRBox()}
    {
        set();
    }

    UsableIonsPopulation(UsableIonsPopulation const& that)
        : Super{pop_dict(that.name())}
        , particleDensity{that.particleDensity}
        , chargeDensity{that.chargeDensity}
        , F{that.F}
        , M{that.M}
        , particles{that.particles}
    {
        set();
    }

    Super& view() { return *this; }
    Super const& view() const { return *this; }
    auto& operator*() { return view(); }
    auto& operator*() const { return view(); }

    _defaults::Grid_t particleDensity, chargeDensity;
    _defaults::UsableVecField_t F;
    _defaults::UsableTensorField_t M;
    UsableParticlesPopulation<ParticleArray_t> particles;
};


template<typename _defaults>
class UsableIons
    : public Ions<typename _defaults::IonPopulation_t, typename _defaults::GridLayout_t>
{
    using GridLayout_t = _defaults::GridLayout_t;
    using Super        = Ions<typename _defaults::IonPopulation_t, GridLayout_t>;

    template<typename PopNames>
    auto static super(PopNames const& pop_names, std::size_t const ppc)
    {
        initializer::PHAREDict dict;
        dict["nbrPopulations"] = pop_names.size();
        for (std::size_t i = 0; i < pop_names.size(); ++i)
            dict["pop" + std::to_string(i)] = pop_dict(pop_names[i], ppc);
        return dict;
    }

    auto static super(Super const supe)
    {
        initializer::PHAREDict dict;
        dict["nbrPopulations"] = supe.size();
        for (std::size_t i = 0; i < supe.size(); ++i)
            dict["pop" + std::to_string(i)] = pop_dict(supe[i].name());
        return dict;
    }

    void set()
    {
        auto&& [_bV, _M, _cd, _md] = Super::getCompileTimeResourcesViewList();
        Vi.set_on(_bV);
        M.set_on(_M);
        _cd.setBuffer(&chargeDensity);
        _md.setBuffer(&massDensity);

        auto& super_pops = Super::getRunTimeResourcesViewList();
        super_pops.clear();
        for (auto& pop : populations)
            super_pops.emplace_back(*pop);
    }

public:
    UsableIons(GridLayout_t const& layout, initializer::PHAREDict const& dict)
        : Super{dict}
        , massDensity{"massDensity", layout, HybridQuantity::Scalar::rho}
        , chargeDensity{"chargeDensity", layout, HybridQuantity::Scalar::rho}
        , Vi{"bulkVel", layout, HybridQuantity::Vector::V}
        , M{"momentumTensor", layout, HybridQuantity::Tensor::M}
    {
        auto& super_pops = Super::getRunTimeResourcesViewList();
        populations.reserve(super_pops.size());
        for (std::size_t i = 0; i < super_pops.size(); ++i)
            populations.emplace_back(dict["pop" + std::to_string(i)], layout);
        set();
    }

    UsableIons(GridLayout_t const& layout, std::vector<std::string> const& pop_names,
               std::size_t const ppc = 0)
        : UsableIons{layout, super(pop_names, ppc)}
    {
    }

    UsableIons(GridLayout_t const& layout, std::string const& pop_name, std::size_t const ppc = 0)
        : UsableIons{layout, std::vector<std::string>{pop_name}, ppc}
    {
    }

    UsableIons(UsableIons&& that)
        : Super(super(*that))
        , massDensity{std::move(that.massDensity)}
        , chargeDensity{std::move(that.chargeDensity)}
        , Vi{std::move(that.Vi)}
        , M{std::move(that.M)}
        , populations{std::move(that.populations)}
    {
        set();
    }

    UsableIons(UsableIons const& that)
        : Super(super(*that))
        , massDensity{that.massDensity}
        , chargeDensity{that.chargeDensity}
        , Vi{that.Vi}
        , M{that.M}
        , populations{that.populations}
    {
        set();
    }

    Super& view() { return *this; }
    Super const& view() const { return *this; }
    auto& operator*() { return view(); }
    auto& operator*() const { return view(); }

    _defaults::Grid_t massDensity, chargeDensity;
    _defaults::UsableVecField_t Vi;
    _defaults::UsableTensorField_t M;
    std::vector<UsableIonsPopulation<_defaults>> populations;
};

template<typename ParticleArray_t, std::size_t interp = 1>
using UsableIonsPopulation_t
    = UsableIonsPopulation<UsableIonsDefaultTypes<ParticleArray_t, interp>>;


template<typename ParticleArray_t, std::size_t interp = 1>
using UsableIons_t = UsableIons<UsableIonsDefaultTypes<ParticleArray_t, interp>>;


} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_ION_POPULATIONS_ION_POPULATION_FIXTURES_HPP */
