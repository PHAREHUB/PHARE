#ifndef PHARE_SOLVER_SOLVER_PPC_MODEL_VIEW_HPP
#define PHARE_SOLVER_SOLVER_PPC_MODEL_VIEW_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "amr/solvers/solver.hpp"
#include <amr/resources_manager/amr_utils.hpp>

namespace PHARE::solver
{
template<typename... Vectors>
void assert_equal_sizes([[maybe_unused]] Vectors const&... vectors)
{
    PHARE_DEBUG_DO( //
        auto tup = std::forward_as_tuple(vectors...);
        assert(core::for_N_all<std::tuple_size_v<decltype(tup)>>(
            [&](auto i) { return std::get<i>(tup).size() == std::get<0>(tup).size(); })); //
    )
}

/*Faraday, Ampere, Ohm Transformers are abstraction that, from the solver viewpoint, act as Faraday,
 * Ampere and Ohm algorithms, but take all patch views and hide the way these are processed, for
 * instance to implement a parallelization decomposition*/
template<typename GridLayout>
class FaradayTransformer
{
    using core_type = PHARE::core::Faraday<GridLayout>;

public:
    template<typename GridLayouts, typename VecFields>
    void operator()(GridLayouts const& layouts, VecFields const& B, VecFields const& E,
                    VecFields& Bnew, double dt)
    {
        assert_equal_sizes(B, E, Bnew);
        for (std::size_t i = 0; i < B.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], faraday_);
            faraday_(*B[i], *E[i], *Bnew[i], dt);
        }
    }

    core_type faraday_;
};

template<typename GridLayout>
class AmpereTransformer
{
    using core_type = PHARE::core::Ampere<GridLayout>;

public:
    template<typename GridLayouts, typename VecFields>
    void operator()(GridLayouts const& layouts, VecFields const& B, VecFields& J)
    {
        assert_equal_sizes(B, J);
        for (std::size_t i = 0; i < B.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], ampere_);
            ampere_(*B[i], *J[i]);
        }
    }
    core_type ampere_;
};

template<typename GridLayout>
class OhmTransformer
{
    using core_type = PHARE::core::Ohm<GridLayout>;

public:
    explicit OhmTransformer(initializer::PHAREDict const& dict)
        : ohm_{dict}
    {
    }

    template<typename GridLayouts, typename VecFields, typename Fields>
    void operator()(GridLayouts const& layouts, Fields const& n, VecFields const& Ve,
                    Fields const& Pe, VecFields const& B, VecFields const& J, VecFields& Enew)
    {
        assert_equal_sizes(n, Ve, Pe, B, J, Enew);
        for (std::size_t i = 0; i < B.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], ohm_);
            ohm_(*n[i], *Ve[i], *Pe[i], *B[i], *J[i], *Enew[i]);
        }
    }

    core_type ohm_;
};



template<typename HybridModel_>
class HybridPPCModelView : public ISolverModelView
{
public:
    using This             = HybridPPCModelView<HybridModel_>;
    using HybridModel_t    = HybridModel_;
    using HybridModel_args = HybridModel_t::type_list::Tuple;
    using IPhysicalModel_t = HybridModel_t::Interface;
    using patch_t          = HybridModel_t::patch_t;
    using Electrons        = HybridModel_t::electrons_t;
    using hierarchy_t      = HybridModel_t::amr_types::hierarchy_t;
    using level_t          = HybridModel_t::amr_types::level_t;
    using Electromag       = HybridModel_t::electromag_type;
    using Ions             = HybridModel_t::ions_type;
    using ParticleArray_t  = Ions::particle_array_type;
    using Particle_t       = ParticleArray_t::value_type;
    using VecFieldT        = HybridModel_t::vecfield_type;
    using FieldT           = HybridModel_t::field_type;
    using GridLayout       = HybridModel_t::gridlayout_type;
    using Faraday_t        = FaradayTransformer<GridLayout>;
    using Ampere_t         = AmpereTransformer<GridLayout>;
    using Ohm_t            = OhmTransformer<GridLayout>;

    struct PatchState_t;

    using value_type = PatchState_t;

    template<bool isConst = false>
    struct iterator;

    HybridPPCModelView(hierarchy_t const& hierarchy, level_t& level, IPhysicalModel_t& model)
        : model_{dynamic_cast<HybridModel_t&>(model)}
    {
        onRegrid(hierarchy, level, model_);
    }

    void onRegrid(hierarchy_t const& hierarchy, level_t& level, HybridModel_t& hybridModel);

    auto begin() { return iterator</*const=*/false>{*this}; }
    auto begin() const { return iterator</*const=*/true>{*this}; }

    auto end() { return iterator</*const=*/false>{*this, states.size()}; }
    auto end() const { return iterator</*const=*/true>{*this, states.size()}; }

    auto& model() { return model_; }
    auto& model() const { return model_; }

    HybridModel_t& model_;
    std::vector<core::aggregate_adapter<PatchState_t>> states;

    /*these vectors allow to access patch data per quantity
       they are kept in synced with "states" so that client code (e.g. Faraday etc.)
       can still have signatures revealing which quantities they use, and not take the states vector
       entirely */
    std::vector<VecFieldT*> electromag_E;
    std::vector<VecFieldT*> electromag_B;
    std::vector<VecFieldT*> electromagPred_E;
    std::vector<VecFieldT*> electromagPred_B;
    std::vector<VecFieldT*> electromagAvg_E;
    std::vector<VecFieldT*> electromagAvg_B;
    std::vector<VecFieldT*> J;
    std::vector<Ions*> ions;
    std::vector<GridLayout*> layouts;
    std::vector<FieldT*> Pe;
    std::vector<FieldT*> N;
    std::vector<VecFieldT*> Ve;

private:
    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};

    auto vecs_tuple()
    {
        return std::forward_as_tuple(electromag_E, electromag_B, electromagPred_E, electromagPred_B,
                                     electromagAvg_E, electromagAvg_B, Pe, N, Ve, J, ions, layouts);
    }
};


template<typename HybridModel>
void HybridPPCModelView<HybridModel>::onRegrid(hierarchy_t const& hierarchy, level_t& level,
                                               HybridModel_t& hybridModel)
{
    auto& hybridState = hybridModel.state;
    auto& rm          = *hybridModel.resourcesManager;

    states.clear();

    for (auto& patch : level)
    {
        {
            auto _ = rm.setOnPatch(*patch, hybridState, electromagPred_, electromagAvg_);
            states.emplace_back(                                             //
                PHARE::amr::layoutFromPatch<GridLayout>(*patch, &hierarchy), //
                hybridState.ions,                                            //
                hybridState.J,                                               //
                hybridState.electromag,                                      //
                electromagPred_,                                             //
                electromagAvg_,                                              //
                hybridState.electrons,                                       //
                patch                                                        //
            );
        }
        assert(states.back().isUsable());
    }
    core::apply(vecs_tuple(), [&](auto& vec) { vec.reserve(states.size()); });
    for (auto& state : states)
    {
        layouts.emplace_back(&state.layout);
        ions.emplace_back(&state.ions);
        J.emplace_back(&state.J);
        Pe.emplace_back(&state.electrons.pressure());
        Ve.emplace_back(&state.electrons.velocity());
        N.emplace_back(&state.electrons.density());
        electromag_E.emplace_back(&state.electromag.E);
        electromag_B.emplace_back(&state.electromag.B);
        electromagAvg_E.emplace_back(&state.electromagAvg.E);
        electromagAvg_B.emplace_back(&state.electromagAvg.B);
        electromagPred_E.emplace_back(&state.electromagPred.E);
        electromagPred_B.emplace_back(&state.electromagPred.B);
    }
}


template<typename HybridModel>
template<bool isConst>
struct HybridPPCModelView<HybridModel>::iterator
{
    using View = std::conditional_t<!isConst, HybridPPCModelView<HybridModel>,
                                    HybridPPCModelView<HybridModel> const>;

    iterator(View& view_, std::size_t const& i = 0)
        : view{view_}
        , idx{i}
    {
    }

    auto& operator*()
    {
        assert(view.states[idx].isUsable());
        return view.states[idx];
    }
    auto& operator*() const
    {
        assert(view.states[idx].isUsable());
        return view.states[idx];
    }

    auto& operator++()
    {
        idx += 1;
        return *this;
    }

    bool operator!=(iterator const& that) const { return &view != &that.view or idx != that.idx; }

    View& view;
    std::size_t idx = 0;
};


template<typename HybridModel>
struct HybridPPCModelView<HybridModel>::PatchState_t
{
    GridLayout layout;
    Ions ions;
    VecFieldT J;
    Electromag electromag, electromagPred, electromagAvg;
    Electrons electrons;
    std::shared_ptr<patch_t> patch;


    bool isUsable() const
    {
        using Tuple               = decltype((*this)());
        auto constexpr tuple_size = std::tuple_size_v<Tuple>;
        Tuple tup                 = (*this)();

        return core::for_N_all<tuple_size>([&](auto i) {
            auto const& [k, v] = std::get<i>(tup);
            auto isUsable      = v->isUsable();
            if (!isUsable)
            {
                PHARE_LOG_LINE_STR("NOT USABLE: " << k);
            }
            return isUsable;
        });
    }

    auto operator()() const
    {
        return core::make_named_tuple(                        //
            std::make_pair("J", &J),                          //
            std::make_pair("ions", &ions),                    //
            std::make_pair("electrons", &electrons),          //
            std::make_pair("electromag", &electromag),        //
            std::make_pair("electromagAvg", &electromagAvg),  //
            std::make_pair("electromagPred", &electromagPred) //
        );
    }
};



} // namespace PHARE::solver



#endif /* PHARE_SOLVER_SOLVER_PPC_MODEL_VIEW_HPP */
