#ifndef PHARE_ION_UPDATER_MULTI_PC_HPP
#define PHARE_ION_UPDATER_MULTI_PC_HPP

#include "ion_updater_def.hpp"

#include "core/numerics/pusher/multi_boris.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/interpolator/interpolating.hpp"



namespace PHARE::core
{


template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdaterMultiPC
{
    using Particles                     = typename Ions::particle_array_type;
    bool constexpr static atomic_interp = Particles::alloc_mode == AllocatorMode::GPU_UNIFIED;

public:
    auto static constexpr dimension    = GridLayout::dimension;
    auto static constexpr interp_order = GridLayout::interp_order;
    using Box_t                        = Box<int, dimension>;
    using Interpolator_t               = Interpolator<dimension, interp_order, atomic_interp>;
    using Interpolating_t              = Interpolating<Particles, interp_order, atomic_interp>;

    using Pusher_t = MultiBorisPusher<GridLayout, Particles, Electromag, Interpolator_t>;

private:
    Interpolating_t interpolator_;

public:
    IonUpdaterMultiPC() {}

    template<typename ModelAccessor_t, typename Boxing_t>
    void updatePopulations(ModelAccessor_t&, Boxing_t const&, double const&,
                           UpdaterMode = UpdaterMode::all);

    static void updateIons(Ions& ions)
    {
        ions.computeChargeDensity();
        ions.computeBulkVelocity();
    }

    void reset() {}


private:
    template<typename ModelAccessor_t, typename Boxing_t>
    void updateAndDepositDomain_(ModelAccessor_t&, Boxing_t const&);

    template<typename ModelAccessor_t, typename Boxing_t>
    void updateAndDepositAll_(ModelAccessor_t&, Boxing_t const&);

    double dt_ = 0;
};




template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelAccessor_t, typename Boxing_t>
void IonUpdaterMultiPC<Ions, Electromag, GridLayout>::updatePopulations(ModelAccessor_t& accessor,
                                                                        Boxing_t const& boxings,
                                                                        double const& dt,
                                                                        UpdaterMode mode)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiPC::updatePopulations");

    for (std::size_t i = 0; i < accessor.size(); ++i)
    {
        auto view      = accessor[i];
        auto [ions, _] = view.args;
        resetMoments(ions);
    }
    dt_ = dt;
    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(accessor, boxings);
    else
        updateAndDepositAll_(accessor, boxings);
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelAccessor_t, typename Boxing_t>
void IonUpdaterMultiPC<Ions, Electromag, GridLayout>::updateAndDepositDomain_(
    ModelAccessor_t& accessor, Boxing_t const& boxings)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiPC::updateAndDepositDomain_");

    if (accessor.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU
    MultiBoris<ModelAccessor_t> in{dt_, accessor};
    Pusher_t::template move<MultiBorisMode::COPY>(in, boxings);
    in.streamer.join();
    in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositDomain_.txt");
#else
    throw std::runtime_error("No available implementation");
#endif
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename ModelAccessor_t, typename Boxing_t>
void IonUpdaterMultiPC<Ions, Electromag, GridLayout>::updateAndDepositAll_(
    ModelAccessor_t& accessor, Boxing_t const& boxings)
{
    PHARE_LOG_SCOPE(1, "IonUpdaterMultiPC::updateAndDepositAll_");

    if (accessor.size() == 0)
        return;

#if PHARE_HAVE_MKN_GPU
    MultiBoris<ModelAccessor_t> in{dt_, accessor};

    Pusher_t::move(in, boxings);

    in.streamer.host([&](auto const i) mutable {
        auto view                 = accessor[i];
        auto [ions, _]            = view.args;
        auto const patch_id       = view.patchID();
        auto const& patch_boxings = boxings.at(patch_id);

        auto const per_pop = [&](auto& pop) {
            auto& domain = pop.domainParticles();
            delete_particles_not_in(domain, patch_boxings.nonLevelGhostBox);
            move_in_ghost_layer(pop.patchGhostParticles(), domain, patch_boxings.domainBox,
                                patch_boxings.nonLevelGhostBox);
            move_in_domain(domain, pop.levelGhostParticles(), patch_boxings.domainBox);
            delete_particles_not_in(pop.levelGhostParticles(), patch_boxings.ghostBox);
            delete_particles_not_in(domain, patch_boxings.domainBox);
        };

        for (auto& pop : ions)
            per_pop(pop);
    });

    in.streamer.host([&](auto const i) mutable {
        auto view            = accessor[i];
        auto [ions, _]       = view.args;
        auto const patch_id  = view.patchID();
        auto const& boxing_i = boxings.at(patch_id);
        Interpolating_t interp;
        for (auto& pop : ions)
        {
            interp.particleToMesh(pop.domainParticles(), boxing_i.layout, pop.particleDensity(),
                                  pop.chargeDensity(), pop.flux());
            interp.particleToMesh(pop.patchGhostParticles(), boxing_i.layout, pop.particleDensity(),
                                  pop.chargeDensity(), pop.flux());
        }
    });

    in.streamer.join();
    in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositAll_.txt");

#else
    throw std::runtime_error("No available implementation");
#endif
}


} // namespace PHARE::core


#endif // PHARE_ION_UPDATER_MULTI_PC_HPP
