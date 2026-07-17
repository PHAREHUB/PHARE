#ifndef PHARE_ION_UPDATER_PP_HPP
#define PHARE_ION_UPDATER_PP_HPP


#include "ion_updater_def.hpp"
#include "core/numerics/pusher/any_boris.hpp"
#include "core/data/particles/particle_array_partitionner.hpp"

#include "core/numerics/interpolator/interpolating.hpp"

namespace PHARE::core
{
template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdaterPP
{
    using Particles                     = Ions::particle_array_type;
    bool constexpr static atomic_interp = Particles::alloc_mode == AllocatorMode::GPU_UNIFIED;
    // static_assert(Particles::alloc_mode == AllocatorMode::GPU_UNIFIED);
    using Particles_vt = Particles::view_t;

public:
    auto static constexpr dimension    = GridLayout::dimension;
    auto static constexpr interp_order = GridLayout::interp_order;
    using Box_t                        = Box<int, dimension>;
    using ParticleArray_t              = Particles;
    using Interpolator_t               = Interpolator<dimension, interp_order, atomic_interp>;
    using Interpolating_t              = Interpolating<Particles_vt, interp_order, atomic_interp>;
    // using Particle_t    = typename Particles::Particle_t;

    using Pusher_t = AnyBorisPusher<dimension, Interpolator_t, GridLayout>;

private:
    Interpolating_t interpolator_;
    // Pusher pusher_;

public:
    IonUpdaterPP() {}

    template<typename Boxing_t>
    void updatePopulations(Ions& ions, Electromag const& em, Boxing_t const& layout,
                           double const& dt, UpdaterMode = UpdaterMode::all);


    static void updateIons(Ions& ions);

    void reset()
    {
        // clear memory
        tmp_particles_ = std::move(make_particles<Particles>(Box_t{}));
    }


private:
    void updateAndDepositDomain_(Ions& ions, Electromag const& em, GridLayout const& layout);
    void updateAndDepositAll_(Ions& ions, Electromag const& em, GridLayout const& layout);


    // dealloced on regridding/load balancing coarsest
    Particles tmp_particles_;

    double dt_ = 0;
};




template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdaterPP<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions, Electromag const& em,
                                                                   Boxing_t const& boxing,
                                                                   double const& dt,
                                                                   UpdaterMode mode)
{
    PHARE_LOG_SCOPE(3, "IonUpdaterPP::updatePopulations");
    resetMoments(ions);
    dt_ = dt;
    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(ions, em, boxing.layout);
    else
        updateAndDepositAll_(ions, em, boxing.layout);
}



template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdaterPP<Ions, Electromag, GridLayout>::updateIons(Ions& ions)
{
    // fixMomentGhosts(ions, layout);
    ions.computeChargeDensity();
    ions.computeBulkVelocity();
}


template<typename Particles, typename GridLayout>
auto pp_partition_particles(Particles& particles, GridLayout const& layout)
{
    return ParticleArrayPartitioner<Particles::alloc_mode, Particles>{particles}(
        std::array{layout.AMRBox(), grow(layout.AMRBox(), GridLayout::nbrParticleGhosts())});
}



template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdaterPP<Ions, Electromag, GridLayout>::updateAndDepositDomain_(Ions& ions,
                                                                         Electromag const& em,
                                                                         GridLayout const& layout)
{
    PHARE_LOG_SCOPE(3, "IonUpdaterPP::updateAndDepositDomain_");

    for (auto& pop : ions)
    {
        auto& domain = (tmp_particles_ = pop.domainParticles());

        Pusher_t pusher{layout, layout.meshSize(), dt_, pop.mass()};
        {
            pusher(domain, em);
            auto iterators = pp_partition_particles(domain, layout);
            domain.resize(iterators[0].size());
            auto v = domain.view();
            interpolator_.particleToMesh(v, layout, pop.density(), pop.flux());
        }

        auto pushGhosts = [&](auto const& inputArray, bool copyInDomain = false) {
            tmp_particles_ = inputArray;
            pusher(tmp_particles_, em);
            auto iterators = pp_partition_particles(tmp_particles_, layout);
            tmp_particles_.resize(iterators[0].size());
            auto v = tmp_particles_.view();
            interpolator_.particleToMesh(v, layout, pop.density(), pop.flux());
        };


        pushGhosts(pop.levelGhostParticles());
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief IonUpdaterPP<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments and particles from time n to n+1
 */
void IonUpdaterPP<Ions, Electromag, GridLayout>::updateAndDepositAll_(Ions& ions,
                                                                      Electromag const& em,
                                                                      GridLayout const& layout)
{
    PHARE_LOG_SCOPE(3, "IonUpdaterPP::updateAndDepositAll_");

    for (auto& pop : ions)
    {
        assert(pop.domainParticles().size());

        auto& domain = pop.domainParticles();
        Pusher_t pusher{layout, layout.meshSize(), dt_, pop.mass()};
        assert(pop.domainParticles().size());
        {
            pusher(domain, em);
            assert(pop.domainParticles().size());
            auto iterators = pp_partition_particles(domain, layout);
            assert(pop.domainParticles().size());
            domain.resize(iterators[0].size());
            assert(pop.domainParticles().size());
        }
        assert(pop.domainParticles().size());

        auto pushGhosts = [&](auto& particles) {
            pusher(particles, em);
            auto iterators = pp_partition_particles(particles, layout);
            domain.reserve(domain.size() + iterators[0].size());
            std::copy(iterators[0].begin(), iterators[0].end(), std::back_inserter(domain));
            particles.erase(iterators.back().end(), particles.end());
        };

        pushGhosts(pop.levelGhostParticles());

        auto v = domain.view();
        interpolator_.particleToMesh(v, layout, pop.density(), pop.flux());
        assert(pop.domainParticles().size());
    }
}



} // namespace PHARE::core


#endif /*PHARE_ION_UPDATER_PP_HPP*/
