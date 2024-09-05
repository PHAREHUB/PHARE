#ifndef _PHARE_AMR_DATA_FIELD_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_
#define _PHARE_AMR_DATA_FIELD_INITIAZILIZERS_SAMRAI_HDF5_INITIALIZER_HPP_

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

#include "amr/data/initializers/samrai_hdf5_initializer.hpp"


namespace PHARE::amr
{


template<typename Field_t, typename GridLayout>
class SamraiHDF5FieldInitializer // no virtual classes needed (yet)
{
public:
    static constexpr auto dimension = GridLayout::dimension;

    SamraiHDF5FieldInitializer() {}

    void load(Field_t& Field, GridLayout const& layout) const;
};



template<typename Field_t, typename GridLayout>
void SamraiHDF5FieldInitializer<Field_t, GridLayout>::load(Field_t& field,
                                                           GridLayout const& layout) const
{
    PHARE_LOG_LINE_STR("SamraiHDF5FieldInitializer::loadParticles");
}



} // namespace PHARE::amr


#endif
