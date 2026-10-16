#ifndef PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP

#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "amr/physical_models/models.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/physical_models/hybrid_model.hpp"


namespace PHARE::solver
{



template<typename Model>
class FaradayLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::Faraday<GridLayout>;

public:
    explicit FaradayLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    void operator()(GridLayout& layout, auto&&... args) { core_type{layout}(args...); }

    void operator()(auto& B, auto& E, auto& Bnew, auto& dt)
    {
        amr::onLevel(amr::make_model_level_accessor(level_, model_, B, E, Bnew), [&](auto& view) {
            auto& [B_v, E_v, Bnew_v] = view.args;
            (*this)(view.layout, B_v, E_v, Bnew_v, dt);
        });
    }

    level_t& level_;
    Model& model_;
};
template<typename Model>
FaradayLevelTransformer(typename Model::amr_types::level_t&, Model&)
    -> FaradayLevelTransformer<Model>;




template<typename Model, core::AmpereBox box = core::AmpereBox{}>
class AmpereLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::Ampere<GridLayout, box>;

public:
    explicit AmpereLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    void operator()(GridLayout& layout, auto&&... args) { core_type{layout}(args...); }

    void operator()(auto& B, auto& J)
    {
        amr::onLevel(amr::make_model_level_accessor(level_, model_, B, J), [&](auto& view) {
            auto& [B_v, J_v] = view.args;
            (*this)(view.layout, B_v, J_v);
        });
    }

    level_t& level_;
    Model& model_;
};


template<typename Model>
AmpereLevelTransformer(typename Model::amr_types::level_t&, Model&)
    -> AmpereLevelTransformer<Model>;




template<typename level_t, typename Model>
struct TimeSetter
{
    void operator()(auto&... quantities)
    {
        amr::onLevel(amr::make_model_level_accessor(level, model, quantities...),
                     [&](auto& view) { view.setTime(newTime, quantities...); });
    }

    level_t& level;
    Model& model;
    double newTime;
};

template<typename level_t, typename Model>
TimeSetter(level_t&, Model&, double) -> TimeSetter<level_t, Model>;



template<typename Model>
struct FieldEvolverDispatchers
{
    // hybrid: J is never communicated — Ampere computes it on physical+1,
    // exactly the layer Ohm reads. MHD keeps ghostbox-1 (its working state).
    static constexpr core::AmpereBox ampere_box
        = is_hybrid_model_v<Model> ? core::AmpereBox{core::AmpereMode::GrownPhysical, 1}
                                   : core::AmpereBox{};

    using Faraday_t = FaradayLevelTransformer<Model>;
    using Ampere_t  = AmpereLevelTransformer<Model, ampere_box>;
};


} // namespace PHARE::solver



#endif /* PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP */
