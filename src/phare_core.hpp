#ifndef PHARE_CORE_INCLUDE_HPP
#define PHARE_CORE_INCLUDE_HPP

#include "core/data/ions/ions.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/electrons/electrons.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/numerics/reconstructions/reconstruction_nghosts.hpp"

#include "phare_simulator_options.hpp"

#include "dict.hpp"

#include <string>
#include <functional>
#include <unordered_map>

namespace PHARE::core
{

template<SimOpts opts>
struct PHARE_Types
{
    auto static constexpr dimension    = opts.dimension;
    auto static constexpr interp_order = opts.interp_order;
    static constexpr auto mhd_reconstruction_nghosts
        = MHDOpts::reconstruction_nghosts_v<opts.reconstruction_type>;

    using Array_t     = NdArrayVector<dimension>;
    using ArrayView_t = NdArrayView<dimension>;

    // Hybrid
    using Grid_t           = Grid<Array_t, HybridQuantity::Scalar>;
    using Field_t          = Field<dimension, HybridQuantity::Scalar>;
    using VecField_t       = VecField<Field_t, HybridQuantity>;
    using SymTensorField_t = SymTensorField<Field_t, HybridQuantity>;
    using Electromag_t     = Electromag<VecField_t>;
    using YeeLayout_t      = GridLayoutImplYee<dimension, interp_order>;
    using GridLayout_t     = GridLayout<YeeLayout_t>;

    using Particle_t      = Particle<dimension>;
    using ParticleAoS_t   = ParticleArray<dimension>;
    using ParticleArray_t = ParticleAoS_t;
    using ParticleSoA_t   = ContiguousParticles<dimension>;

    using MaxwellianParticleInitializer_t
        = MaxwellianParticleInitializer<ParticleArray_t, GridLayout_t>;
    using IonPopulation_t = IonPopulation<ParticleArray_t, VecField_t, SymTensorField_t>;
    using Ions_t          = Ions<IonPopulation_t, GridLayout_t>;
    using Electrons_t     = Electrons<Ions_t>;

    using ParticleInitializerFactory_t = ParticleInitializerFactory<ParticleArray_t, GridLayout_t>;


    using Grid_MHD     = Grid<Array_t, MHDQuantity::Scalar>;
    using Field_MHD    = Field<dimension, MHDQuantity::Scalar>;
    using VecField_MHD = VecField<Field_MHD, MHDQuantity>;

    using YeeLayout_MHD = GridLayoutImplYeeMHD<dimension, interp_order, mhd_reconstruction_nghosts>;
    using GridLayout_MHD = GridLayout<YeeLayout_MHD>;
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
