#ifndef _PHARE_CORE_DATA_IONS_PARTICLE_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_CORE_DATA_IONS_PARTICLE_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_

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

#include "hdf5/detail/h5/h5_file.hpp"


#include "SAMRAI/hier/PatchDataRestartManager.h"


namespace PHARE::amr
{


template<typename ParticleArray, typename GridLayout>
class SamraiH5Interface
{
public:
    static SamraiH5Interface& INSTANCE()
    {
        static SamraiH5Interface i;
        return i;
    }

    void populate_from(std::string const& dir, int const& idx, int const& mpi_size);

    NO_DISCARD auto static getRestartFileFullPath(std::string path, int const& idx,
                                                  int const& mpi_size, int const& rank)
    {
        return path                                                          //
               + "/restore." + SAMRAI::tbox::Utilities::intToString(idx, 6)  //
               + "/nodes." + SAMRAI::tbox::Utilities::nodeToString(mpi_size) //
               + "/proc." + SAMRAI::tbox::Utilities::processorToString(rank);
    }


private:
    std::unordered_map<std::string, std::string> box2dataset;
};


/*
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/d_box
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/d_ghost_box
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/d_ghosts
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/d_timestamp
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/domainParticles_charge
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/domainParticles_delta
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/domainParticles_iCell
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/domainParticles_v
/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/protons##default/domainParticles_weight
*/


template<typename ParticleArray, typename GridLayout>
void SamraiH5Interface<ParticleArray, GridLayout>::populate_from(std::string const& dir,
                                                                 int const& idx,
                                                                 int const& mpi_size)
{
    for (int rank = 0; rank < mpi_size; ++rank)
    {
        auto const hdf5_filepath = getRestartFileFullPath(dir, idx, mpi_size, rank);

        hdf5::h5::HighFiveFile h5File{hdf5_filepath, HighFive::File::ReadOnly, /*para=*/false};

        PHARE_LOG_LINE_STR("SamraiH5Interface::populate_from");

        auto groups = h5File.scan_for_groups({"level_0000", "domainParticles_charge"});

        for (auto const& g : groups)
        {
            PHARE_LOG_LINE_STR(g);
        }
    }
}



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
