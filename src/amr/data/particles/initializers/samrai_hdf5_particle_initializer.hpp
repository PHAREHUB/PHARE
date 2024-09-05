#ifndef _PHARE_AMR_DATA_PARTICLE_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_PARTICLE_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_

#include <memory>
#include <random>
#include <cassert>
#include <functional>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ions/particle_initializers/particle_initializer.hpp"
#include "core/data/particles/particle.hpp"
#include "initializer/data_provider.hpp"
#include "core/utilities/point/point.hpp"
#include "core/def.hpp"
#include "core/logger.hpp"

#include "hdf5/detail/h5/h5_file.hpp"


#include "SAMRAI/hier/PatchDataRestartManager.h"


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
    PHARE_LOG_LINE_STR("SamraiHDF5ParticleInitializer::loadParticles");
}




} // namespace PHARE::amr


#endif
