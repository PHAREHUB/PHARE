#ifndef PHARE_BENCH_CORE_BENCH_H
#define PHARE_BENCH_CORE_BENCH_H

#include "mkn/kul/dbg.hpp"
#include "mkn/kul/log.hpp"

#include "benchmark/benchmark.h"

#include "phare_core.h"
#include "core/utilities/types.h"

namespace PHARE::core::bench
{
template<std::size_t dim>
using Field = PHARE::core::Field<PHARE::core::NdArrayVector<dim>,
                                 typename PHARE::core::HybridQuantity::Scalar>;
template<std::size_t dim>
using VecField
    = PHARE::core::VecField<PHARE::core::NdArrayVector<dim>, typename PHARE::core::HybridQuantity>;

// clang-format off
template<std::size_t dim>
PHARE::core::Particle<dim> particle(int icell = 15)
{
    return {
        /*.weight = */ .00001,
        /*.charge = */ .00001,
        /*.iCell  = */ PHARE::core::ConstArray<int, dim>(icell),
        /*.delta  = */ PHARE::core::ConstArray<double, dim>(.5),
        /*.v      = */ {{.00001, .00001, .00001}},
        /* E xyz      = */ .00001, .00001, .00001,
        /* B xyz      = */ .00001, .00001, .00001,
    };
}
// clang-format on

template<std::size_t dim>
PHARE::core::ParticleArray<dim> make_particles(std::size_t ppc)
{
    return {ppc, particle<dim>()};
}

template<typename Particles, typename Point>
void disperse(Particles& particles, Point lo, Point up, std::optional<int> seed = std::nullopt)
{
    auto gen = [&]() {
        if (!seed.has_value())
        {
            std::random_device rd;
            std::seed_seq seed_seq{rd(), rd(), rd(), rd(), rd(), rd(), rd()};
            return std::mt19937_64(seed_seq);
        }
        return std::mt19937_64(*seed);
    }();
    for (std::size_t i = 0; i < Particles::dimension; i++)
    {
        std::uniform_int_distribution<> distrib(lo[i], up[i]);
        for (auto& particle : particles)
            particle.iCell[i] = distrib(gen);
    }
}
template<typename Particles>
void disperse(Particles& particles, std::size_t lo, std::size_t up,
              std::optional<int> seed = std::nullopt)
{
    auto static dim = Particles::dimension;

    disperce(particles, core::ConstArray<int, dim>(lo), core::ConstArray<int, dim>(up), seed);
}

template<std::size_t dim, typename Box>
auto make_particles(std::size_t ppc, Box disperse_in, std::optional<int> seed = std::nullopt)
{
    auto particles = make_particles<dim>(ppc);
    disperse(particles, disperse_in.lower, disperse_in.upper - 1u, seed);
    return particles;
}

template<typename GridLayout, typename Quantity, std::size_t dim = GridLayout::dimension>
Field<dim> field(std::string key, Quantity type, GridLayout const& layout)
{
    Field<dim> feeld{key, type, layout.allocSize(type)};
    std::fill(feeld.begin(), feeld.end(), 1);
    return feeld;
}


template<typename GridLayout, std::size_t dim = GridLayout::dimension>
auto rho(GridLayout const& layout)
{
    return field("rho", HybridQuantity::Scalar::rho, layout);
}


template<typename GridLayout>
class _VF_ : public VecField<GridLayout::dimension>
{
public:
    using Super = VecField<GridLayout::dimension>;

    _VF_(_VF_ const&) = delete;
    _VF_(_VF_&&)      = delete;
    auto& operator=(_VF_ const&) = delete;
    auto& operator=(_VF_&&) = delete;

    _VF_(GridLayout const& layout, HybridQuantity::Vector v_id, std::string id,
         std::array<HybridQuantity::Scalar, 3> quantities)
        : Super{id, v_id}
        , xyz{field(id + "x", quantities[0], layout), field(id + "y", quantities[1], layout),
              field(id + "z", quantities[2], layout)}
    {
        Super::setBuffer(id + "_x", &xyz[0]);
        Super::setBuffer(id + "_y", &xyz[1]);
        Super::setBuffer(id + "_z", &xyz[2]);
    }


    _VF_(GridLayout const& layout, HybridQuantity::Vector v_id, std::string id)
        : _VF_{layout,
               v_id,
               id,
               {HybridQuantity::Scalar::Vx, HybridQuantity::Scalar::Vy, HybridQuantity::Scalar::Vz}}
    {
    }

private:
    std::array<Field<GridLayout::dimension>, 3> xyz;
};

template<typename GridLayout>
class Flux : public _VF_<GridLayout>
{
public:
    using Super = _VF_<GridLayout>;

    Flux(GridLayout const& layout)
        : Super{layout, HybridQuantity::Vector::V, "F"}
    {
    }
};

template<typename GridLayout>
class J : public _VF_<GridLayout>
{
public:
    using Super = _VF_<GridLayout>;

    J(GridLayout const& layout)
        : Super{layout, HybridQuantity::Vector::J, "J"}
    {
    }
};


template<typename GridLayout>
struct Electromag
{
    Electromag(Electromag const&) = delete;

    Electromag(GridLayout const& layout)
        : E{layout,
            HybridQuantity::Vector::E,
            "E",
            {HybridQuantity::Scalar::Ex, HybridQuantity::Scalar::Ey, HybridQuantity::Scalar::Ez}}
        , B{layout,
            HybridQuantity::Vector::B,
            "B",
            {HybridQuantity::Scalar::Bx, HybridQuantity::Scalar::By, HybridQuantity::Scalar::Bz}}
    {
    }

    auto operator()() { return std::forward_as_tuple(E, B); }

    _VF_<GridLayout> E, B;
};


auto inline ions_dict(std::size_t idx = 0, std::string pop = "protons")
{
    PHARE::initializer::PHAREDict ions;
    ions["nbrPopulations"] = 1;
    ions["pop0"]["mass"]   = 1.;
    ions["pop0"]["name"]   = pop + std::to_string(idx);
    // not used, but needed in construction
    ions["pop0"]["particle_initializer"] = std::string("None");
    return ions;
}

template<typename GridLayout_>
struct HybridPatch
{
    using GridLayout = GridLayout_;

    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using PHARE_TYPES     = core::PHARE_Types<dimension, interp_order>;
    using Array_t         = typename PHARE_TYPES::Array_t;
    using Field_t         = typename PHARE_TYPES::Field_t;
    using ParticleArray_t = typename PHARE_TYPES::ParticleArray_t;
    using Ions_t          = typename PHARE_TYPES::Ions_t;
    using ParticlesPack_t = PHARE::core::ParticlesPack<ParticleArray_t>;

    struct Layout
    {
        GridLayout layout;
    };

    struct Particles
    {
        Particles(GridLayout const& layout, std::size_t ppc = 100)
            : domain{core::bench::make_particles<dimension>(ppc, layout.AMRBox())}
            , patch_ghost{0}
            , level_ghost{0}
            , particles_pack{&domain, &patch_ghost, &level_ghost}
        {
        }

        ParticleArray_t domain, patch_ghost, level_ghost;
        ParticlesPack_t particles_pack;
    };


    template<bool view = 0>
    struct Vars
    {
        template<typename T>
        using Wrap_t = std::conditional_t<view, T*, T>;

        template<bool V = view, std::enable_if_t<!V, bool> = 0>
        Vars(GridLayout const& layout, Particles& particles)
            : EM{layout}
            , F{layout}
            , J{layout}
            , rho{core::bench::rho(layout)}
            , ions{ions_dict()}
        {
        }

        template<typename HybridPatchState, bool V = view, std::enable_if_t<V, bool> = 0>
        Vars(HybridPatchState& state)
            : EM{&state.EM}
            , F{&state.F}
            , J{&state.J}
            , rho{&state.rho}
            , ions{&state.ions}
        {
        }

        Wrap_t<core::bench::Electromag<GridLayout>> EM;
        Wrap_t<core::bench::Flux<GridLayout>> F;
        Wrap_t<core::bench::J<GridLayout>> J;
        Wrap_t<Field_t> rho;
        Wrap_t<Ions_t> ions;
    };

    struct View : public Layout, public Vars<1>
    {
        static constexpr auto dimension = GridLayout::dimension;
        using HybridQuantity            = core::HybridQuantity;
        using LocalBox                  = core::Box<std::uint32_t, dimension>;
        using Vars_t                    = Vars<1>;
        using Layout::layout;

        template<typename HybridPatchState>
        View(HybridPatchState& state)
            : Layout{state.layout}
            , Vars_t{state}
            , operate_box{layout.AMRToLocal(layout.AMRBox())}
        {
        }

        auto B_boxes() const { return _boxes(HybridQuantity::B()); }
        auto E_boxes() const { return _boxes(HybridQuantity::E()); }
        auto J_boxes() const { return _boxes(HybridQuantity::J()); }
        auto V_boxes() const { return _boxes(HybridQuantity::V()); }
        auto F_boxes() const { return V_boxes(); }

        auto _boxes(std::array<HybridQuantity::Scalar, 3> quantities) const
        {
            return core::generate([&](auto qty) { return make_box(qty); }, quantities);
        }

        auto make_box(HybridQuantity::Scalar quantity) const
        {
            auto const primal_directions = layout.primal_directions(quantity);
            auto upper                   = operate_box.upper;
            for (std::size_t i = 0; i < dimension; ++i)
                if (upper[i] == local_box.upper[i] and primal_directions[i])
                    upper[i] += 1;
            return LocalBox{operate_box.lower, upper};
        }


        LocalBox const operate_box;
        LocalBox const local_box{layout.AMRToLocal(layout.AMRBox())};
        double const mass            = 1;
        std::size_t const part_start = 0;
        std::size_t const part_end   = 0;
    };



    struct State : public Layout, public Particles, public Vars<0>
    {
        static constexpr auto dimension = GridLayout::dimension;
        using view_t                    = HybridPatch<GridLayout>::View;

        using Vars_t = Vars<0>;
        using Layout::layout;
        using Particles::particles_pack;
        using Vars_t::rho;


        State(State const&) = delete;
        State(State&&)      = delete;
        State(GridLayout& layout_, std::size_t ppc = 100)
            : Layout{layout_}
            , Particles{layout_, ppc}
            , Vars_t{layout_, static_cast<Particles&>(*this)}
            , masses{1}
        {
            auto& pop          = Vars_t::ions.getRunTimeResourcesUserList()[0];
            std::string pop_id = "protons0";
            pop.setBuffer(pop_id, &particles_pack);
            pop.setBuffer(pop_id + "_rho", &rho);

            auto& flux = pop.flux();
            auto& F    = Vars_t::F;
            flux.setBuffer(pop_id + "_flux_x", &F[0]);
            flux.setBuffer(pop_id + "_flux_y", &F[1]);
            flux.setBuffer(pop_id + "_flux_z", &F[2]);
        }

        template<typename Vec>
        static auto make_unique(Vec lo, Vec up, std::size_t ppc = 100)
        {
            assert(lo.size() == up.size() and lo.size() == dimension);
            std::array<std::uint32_t, dimension> cells;
            std::array<double, dimension> meshSize, origin;
            for (std::size_t i = 0; i < cells.size(); ++i)
            {
                assert(up[i] >= lo[i]);
                cells[i]    = up[i] - lo[i] + 1;
                meshSize[i] = 1. / cells[i];
                origin[i]   = lo[i] * meshSize[i];
            }
            GridLayout layout{meshSize, cells, core::Point<double, dimension>{origin}};
            return std::make_unique<State>(layout, ppc);
        }

        std::vector<double> masses;
    };
};




} // namespace PHARE::core::bench

#endif /*PHARE_BENCH_CORE_BENCH_H*/
