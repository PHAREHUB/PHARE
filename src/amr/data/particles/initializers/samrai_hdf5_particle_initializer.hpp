#ifndef _PHARE_AMR_DATA_PARTICLE_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_PARTICLE_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_

#include <memory>
#include <random>
#include <cassert>
#include <functional>

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ions/particle_initializers/particle_initializer.hpp"
#include "core/data/particles/particle.hpp"
#include "initializer/data_provider.hpp"
#include "core/utilities/point/point.hpp"


#include "hdf5/detail/h5/h5_file.hpp"

#include "core/data/particles/particle_packer.hpp"
#include "amr/data/field/initializers/samrai_hdf5_field_initializer.hpp"


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

    auto const& overlaps
        = SamraiH5Interface<GridLayout>::INSTANCE().box_intersections(layout.AMRBox());
    for (auto const& [overlap_box, h5FilePtr, pdataptr] : overlaps)
    {
        auto& h5File              = *h5FilePtr;
        auto& pdata               = *pdataptr;
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
            if (auto const p = soa.copy(i); core::isIn(core::Point{p.iCell}, overlap_box))
                particles.push_back(p);
    }
}




} // namespace PHARE::amr


#endif
