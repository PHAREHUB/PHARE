#ifndef PHARE_CORE_INCLUDE_HPP
#define PHARE_CORE_INCLUDE_HPP

#include "core/def/phare_config.hpp"

#include "core/data/grid/grid.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/electrons/electrons.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/numerics/reconstructions/reconstruction_nghosts.hpp"
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"

#include "phare_simulator_options.hpp"

#include "dict.hpp"

#include <string>
#include <cstddef>
#include <functional>
#include <string_view>
#include <unordered_map>

namespace PHARE::core
{



using field_value_type = double;


template<typename GridLayout_t, auto L_, auto A_>
struct UsingResolver
{
    bool static constexpr c_ordering   = true;
    auto static constexpr dimension    = GridLayout_t::dimension;
    auto static constexpr interp_order = GridLayout_t::interp_order;

    using Array_t = NdArrayVector<dimension, field_value_type, c_ordering, A_>;
    using Grid_t  = Grid<Array_t, HybridQuantity::Scalar>;
    using Field_t = Field<dimension, HybridQuantity::Scalar, field_value_type, A_>;
};

template<typename GridLayout_t, auto A_>
struct UsingResolver<GridLayout_t, LayoutMode::AoSTS, A_>
{
    bool static constexpr c_ordering   = true;
    auto static constexpr dimension    = GridLayout_t::dimension;
    auto static constexpr interp_order = GridLayout_t::interp_order;

    using base_types     = UsingResolver<GridLayout_t, LayoutMode::AoS, A_>;
    using Grid_base_type = base_types::Grid_t;
    using Field_base_type
        = basic::Field<FieldOpts<HybridQuantity::Scalar, field_value_type>{dimension, A_}>;
    using Grid_t  = GridTileSet<GridLayout_t, Grid_base_type, Field_base_type>;
    using Field_t = FieldTileSet<GridLayout_t, Grid_base_type, Field_base_type>;
};

template<typename GridLayout_t, auto A_>
struct UsingResolver<GridLayout_t, LayoutMode::AoSCMTS, A_>
    : UsingResolver<GridLayout_t, LayoutMode::AoSTS, A_>
{
};

template<typename GridLayout_t, auto A_>
struct UsingResolver<GridLayout_t, LayoutMode::AoSPCTS, A_>
    : UsingResolver<GridLayout_t, LayoutMode::AoSTS, A_>
{
};




template<SimOpts opts>
struct PHARE_Types
{
    auto static constexpr dimension    = opts.dimension;
    auto static constexpr interp_order = opts.interp_order;
    static constexpr auto mhd_reconstruction_nghosts
        = MHDOpts::reconstruction_nghosts_v<opts.reconstruction_type>;
    auto static constexpr layout_mode    = opts.layout_mode;
    auto static constexpr allocator_mode = opts.alloc_mode;

    using ParticleSoA_t   = SoAParticleArray<dimension>;
    using ParticleArray_t = ParticleArray<ParticleArrayOptions{
        dimension, layout_mode, StorageMode::VECTOR, allocator_mode}>;

    using YeeLayout_t  = GridLayoutImplYee<dimension, interp_order>;
    using GridLayout_t = GridLayout<YeeLayout_t>;

    using Resolver_t = UsingResolver<GridLayout_t, layout_mode, allocator_mode>;

    using Array_t
        = NdArrayVector<dimension, field_value_type, Resolver_t::c_ordering, allocator_mode>;

    using Grid_t  = Resolver_t::Grid_t;  // layout dependent
    using Field_t = Resolver_t::Field_t; // layout dependent

    using VecField_t       = VecField<Field_t, HybridQuantity>;
    using Electromag_t     = Electromag<VecField_t>;
    using SymTensorField_t = SymTensorField<Field_t, HybridQuantity>;
    using IonPopulation_t  = IonPopulation<ParticleArray_t, VecField_t, SymTensorField_t>;
    using Ions_t           = Ions<IonPopulation_t, GridLayout_t>;
    using Electrons_t      = Electrons<Ions_t>;
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


namespace PHARE::strings
{

constexpr static std::string_view cma = ",";
}

#endif // PHARE_CORE_INCLUDE_HPP
