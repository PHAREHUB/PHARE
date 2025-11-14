#ifndef _PHARE_AMR_DATA_PARTICLE_INITIALIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_PARTICLE_INITIALIZERS_SAMRAI_HDF5_INITIALIZER_HPP_


#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/ions/particle_initializers/particle_initializer.hpp"

#include "amr/data/initializers/samrai_hdf5_initializer.hpp"

#include "core/data/particles/particle_packer.hpp"

#if !PHARE_HAS_HIGHFIVE
#error // HIGHFIVE REQUIRED!
#endif // PHARE_HAS_HIGHFIVE

#include <cassert>

namespace PHARE::amr
{


template<typename ParticleArray, typename GridLayout>
class SamraiHDF5ParticleInitializer : public core::ParticleInitializer<ParticleArray, GridLayout>
{
public:
    static constexpr auto dimension = GridLayout::dimension;



    SamraiHDF5ParticleInitializer() {}



    void loadParticles(ParticleArray& particles, GridLayout const& layout,
                       std::string const& popname) const override;
};



template<typename ParticleArray, typename GridLayout>
void SamraiHDF5ParticleInitializer<ParticleArray, GridLayout>::loadParticles(
    ParticleArray& particles, GridLayout const& layout, std::string const& popname) const
{
    using Packer = core::ParticlePacker<ParticleArray::dimension>;

    auto const& dest_box = layout.AMRBox();
    auto const& overlaps = SamraiH5Interface<GridLayout>::INSTANCE().box_intersections(dest_box);

    for (auto const& [overlap_box, h5FilePtr, pdataptr] : overlaps)
    {
        auto& h5File = *h5FilePtr;
        auto& pdata  = *pdataptr;

        std::string const poppath = pdata.base_path + "/" + popname + "##default/domainParticles_";
        core::ContiguousParticles<ParticleArray::dimension> soa{0};

        {
            std::size_t part_idx = 0;
            core::apply(soa.as_tuple(), [&](auto& arg) {
                auto const datapath = poppath + Packer::keys()[part_idx++];
                h5File.file().getDataSet(datapath).read(arg);
            });
        }

        for (std::size_t i = 0; i < soa.size(); ++i)
            if (auto const p = soa.copy(i); core::isIn(core::Point{p.iCell}, dest_box))
                particles.push_back(p);
    }
}


} // namespace PHARE::amr


#endif /*_PHARE_AMR_DATA_PARTICLE_INITIALIZERS_SAMRAI_HDF5_INITIALIZER_HPP_*/
