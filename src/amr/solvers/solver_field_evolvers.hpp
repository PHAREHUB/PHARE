#ifndef PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP

#include "core/data/field/field_tiles.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "amr/resources_manager/amr_utils.hpp"

namespace PHARE::solver
{

class FaradaySingleTransformer
{
    template<typename GridLayout>
    void operate(GridLayout const& layout, auto&&... args)
    {
        core::Faraday<GridLayout>{layout}(args...);
    }

public:
    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField const& E, VecField& Bnew,
                    double dt)
    {
        using field_type = VecField::field_type;

        if constexpr (core::is_field_tile_set_v<field_type>)
        {
            core::tile_exec_with_layout(
                [&](auto& layout, auto&&... args) { operate(layout, args...); }, B, E, Bnew, dt);

            core::sync_inner_ghosts(Bnew);
        }
        else
        {
            operate(layout, B, E, Bnew, dt);
        }
    }
};

template<typename Model>
class FaradayLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;

public:
    explicit FaradayLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    template<typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField const& E, VecField& Bnew,
                    double dt)
    {
        FaradaySingleTransformer{}(layout, B, E, Bnew, dt);
    }

    void operator()(auto& B, auto& E, auto& Bnew, auto& dt)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, B, E, Bnew))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            (*this)(layout, B, E, Bnew, dt);
        }
    }

    level_t& level_;
    Model& model_;
};

template<typename Model>
FaradayLevelTransformer(typename Model::amr_types::level_t&, Model&)
    -> FaradayLevelTransformer<Model>;


class AmpereSingleTransformer
{
    template<typename GridLayout>
    void operate(GridLayout const& layout, auto&&... args)
    {
        core::Ampere<GridLayout>{layout}(args...);
    }

public:
    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField& J)
    {
        using field_type = VecField::field_type;

        if constexpr (core::is_field_tile_set_v<field_type>)
        {
            core::tile_exec_with_layout(
                [&](auto& layout, auto&&... args) { operate(layout, args...); }, B, J);

            core::sync_inner_ghosts(J);
        }
        else
        {
            operate(layout, B, J);
        }
    }
};

template<typename Model>
class AmpereLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;

public:
    explicit AmpereLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    template<typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField& J)
    {
        AmpereSingleTransformer{}(layout, B, J);
    }

    void operator()(auto& B, auto& J)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, B, J))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            (*this)(layout, B, J);
        }
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
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, quantities...))
            (model.resourcesManager->setTime(quantities, *patch, newTime), ...);
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
    using Faraday_t = FaradayLevelTransformer<Model>;
    using Ampere_t  = AmpereLevelTransformer<Model>;
};

} // namespace PHARE::solver

#endif /* PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP */
