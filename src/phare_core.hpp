
#ifndef PHARE_CORE_INCLUDE_HPP
#define PHARE_CORE_INCLUDE_HPP

#include "core/logger.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/electrons/electrons.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/models/physical_state.hpp"
#include "core/models/physical_state.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/utilities/algorithm.hpp"

#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <unordered_map>

#include "cppdict/include/dict.hpp"

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
    using ParticleSoA_t   = PHARE::core::ParticleArray<dimension, true>;
    using ParticleArray_t = ParticleSoA_t;
    // using ParticleArray_t = ParticleAoS_t;

    using MaxwellianParticleInitializer_t
        = PHARE::core::MaxwellianParticleInitializer<ParticleArray_t, GridLayout_t>;
    using IonPopulation_t = PHARE::core::IonPopulation<ParticleArray_t, VecField_t, GridLayout_t>;
    using Ions_t          = PHARE::core::Ions<IonPopulation_t, GridLayout_t>;
    using Electrons_t     = PHARE::core::Electrons<Ions_t>;

    using ParticleInitializerFactory
        = PHARE::core::ParticleInitializerFactory<ParticleArray_t, GridLayout_t>;
};

struct PHARE_Sim_Types
{
    using SimFunctorParams = cppdict::Dict<int, unsigned int, double, std::size_t>;
    using SimFunctor       = std::function<void(SimFunctorParams const& /*params*/)>;
    using SimulationFunctors // code place id -> function_id -> function
        = std::unordered_map<std::string, std::unordered_map<std::string, SimFunctor>>;
};

} // namespace PHARE::core

#endif // PHARE_CORE_INCLUDE_HPP
