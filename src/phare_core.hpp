
#ifndef PHARE_CORE_INCLUDE_HPP
#define PHARE_CORE_INCLUDE_HPP

#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <unordered_map>

#include "core/logger.hpp"
#include "core/def/types.hpp"

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

#include "cppdict/include/dict.hpp"

namespace PHARE::core
{
template<std::size_t dimension, std::size_t interp_order = 0>
struct CPU_Types;


template<std::size_t dimension>
struct CPU_Types<dimension, 0>
{
    using Array_t         = PHARE::core::NdArrayVector<dimension>;
    using VecField_t      = PHARE::core::VecField<Array_t, PHARE::core::HybridQuantity>;
    using Field_t         = PHARE::core::Field<Array_t, PHARE::core::HybridQuantity::Scalar>;
    using Electromag_t    = PHARE::core::Electromag<VecField_t>;
    using Particle_t      = PHARE::core::Particle<dimension>;
    using ParticleArray_t = PHARE::core::ParticleArray<Particle_t>;
};

template<std::size_t dimension, std::size_t interp_order>
struct CPU_Types
{
    using Array_t         = typename CPU_Types<dimension>::Array_t;
    using VecField_t      = typename CPU_Types<dimension>::VecField_t;
    using Field_t         = typename CPU_Types<dimension>::Field_t;
    using Electromag_t    = typename CPU_Types<dimension>::Electromag_t;
    using Particle_t      = typename CPU_Types<dimension>::Particle_t;
    using ParticleArray_t = typename CPU_Types<dimension>::ParticleArray_t;

    using YeeLayout_t     = PHARE::core::GridLayoutImplYee<dimension, interp_order>;
    using GridLayout_t    = PHARE::core::GridLayout<YeeLayout_t>;
    using IonPopulation_t = PHARE::core::IonPopulation<ParticleArray_t, VecField_t, GridLayout_t>;
    using Ions_t          = PHARE::core::Ions<IonPopulation_t, GridLayout_t>;
    using Electrons_t     = PHARE::core::Electrons<Ions_t>;

    using ParticleInitializerFactory
        = PHARE::core::ParticleInitializerFactory<ParticleArray_t, GridLayout_t>;
    using MaxwellianParticleInitializer_t
        = PHARE::core::MaxwellianParticleInitializer<ParticleArray_t, GridLayout_t>;
};



template<std::size_t dimension, std::size_t interp_order = 0>
struct PHARE_Types;


template<std::size_t dimension_>
struct PHARE_Types<dimension_, 0>
{
    static auto constexpr dimension    = dimension_;
    static auto constexpr interp_order = 0;

    template<typename Type>
    using Allocator_t     = typename PHARE::Vector<Type>::allocator_type;
    using Particle_t      = PHARE::core::Particle<dimension>;
    using ParticleArray_t = PHARE::core::ParticleArray<Particle_t, Allocator_t<Particle_t>>;

    using Array_t
        = PHARE::core::NdArrayVector<dimension, double, /*c_ordering=*/true, Allocator_t<double>>;

    using VecField_t   = PHARE::core::VecField<Array_t, PHARE::core::HybridQuantity>;
    using Field_t      = PHARE::core::Field<Array_t, PHARE::core::HybridQuantity::Scalar>;
    using Electromag_t = PHARE::core::Electromag<VecField_t>;
};

template<std::size_t dimension_, std::size_t interp_order_>
struct PHARE_Types
{
    static auto constexpr dimension    = dimension_;
    static auto constexpr interp_order = interp_order_;

    template<typename Type>
    using Allocator_t     = typename PHARE::Vector<Type>::allocator_type;
    using Particle_t      = PHARE::core::Particle<dimension>;
    using ParticleArray_t = PHARE::core::ParticleArray<Particle_t, Allocator_t<Particle_t>>;
    using Array_t
        = PHARE::core::NdArrayVector<dimension, double, /*c_ordering=*/true, Allocator_t<double>>;
    using VecField_t   = PHARE::core::VecField<Array_t, PHARE::core::HybridQuantity>;
    using Field_t      = PHARE::core::Field<Array_t, PHARE::core::HybridQuantity::Scalar>;
    using Electromag_t = PHARE::core::Electromag<VecField_t>;

    using YeeLayout_t  = PHARE::core::GridLayoutImplYee<dimension, interp_order>;
    using GridLayout_t = PHARE::core::GridLayout<YeeLayout_t>;
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
