
#ifndef PHARE_CORE_INCLUDE_H
#define PHARE_CORE_INCLUDE_H

#include "core/data/electromag/electromag.h"
#include "core/data/electrons/electrons.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/ions/ion_population/ion_population.h"
#include "core/data/ions/ions.h"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "core/data/ndarray/ndarray_vector.h"
#include "core/data/particles/particle_array.h"
#include "core/data/vecfield/vecfield.h"
#include "core/models/physical_state.h"
#include "core/models/physical_state.h"
#include "core/utilities/meta/meta_utilities.h"
#include "core/utilities/algorithm.h"

namespace PHARE::core
{
template<std::size_t dimension_, std::size_t interp_order_>
struct PHARE_Types
{
    static auto constexpr dimension    = dimension_;
    static auto constexpr interp_order = interp_order_;

    using Array_t      = PHARE::core::NdArrayVector<dimension>;
    using VecField_t   = PHARE::core::VecField<Array_t, PHARE::core::HybridQuantity>;
    using Field_t      = PHARE::core::Field<Array_t, PHARE::core::HybridQuantity::Scalar>;
    using Electromag_t = PHARE::core::Electromag<VecField_t>;
    using YeeLayout_t  = PHARE::core::GridLayoutImplYee<dimension, interp_order>;
    using GridLayout_t = PHARE::core::GridLayout<YeeLayout_t>;

    using Particle_t      = PHARE::core::Particle<dimension>;
    using ParticleAoS_t   = PHARE::core::ParticleArray<dimension>;
    using ParticleArray_t = ParticleAoS_t;
    using ParticleSoA_t   = PHARE::core::ContiguousParticles<dimension>;


    using MaxwellianParticleInitializer_t
        = PHARE::core::MaxwellianParticleInitializer<ParticleArray_t, GridLayout_t>;
    using IonPopulation_t = PHARE::core::IonPopulation<ParticleArray_t, VecField_t, GridLayout_t>;
    using Ions_t          = PHARE::core::Ions<IonPopulation_t, GridLayout_t>;
    using Electrons_t     = PHARE::core::Electrons<Ions_t>;

    using ParticleInitializerFactory
        = PHARE::core::ParticleInitializerFactory<ParticleArray_t, GridLayout_t>;
};

} // namespace PHARE::core

#endif // PHARE_CORE_INCLUDE_H
