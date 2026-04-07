#ifndef PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "amr/resources_manager/amr_utils.hpp"


namespace PHARE::solver
{


template<typename Model>
class FaradayTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::Faraday<GridLayout>;

public:
    void operator()(GridLayout& layout, auto&&... args) { core_type{layout}(args...); }

    void operator()(level_t& level, auto& model, auto& B, auto& E, auto& Bnew, auto& dt)
    {
        auto& rm    = *model.resourcesManager;
        auto& state = model.state;
        for (auto& patch : rm.enumerate(level, state, B, E, Bnew))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            (*this)(layout, B, E, Bnew, dt);
        }
    }
};


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

template<typename Model>
class AmpereTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::Ampere<GridLayout>;

public:
    void operator()(GridLayout& layout, auto&&... args) { core_type{layout}(args...); }

    void operator()(level_t& level, auto& model, auto& B, auto& J)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, B, J))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            (*this)(layout, B, J);
        }
    }
};




template<typename Model>
class AmpereLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::Ampere<GridLayout>;

public:
    explicit AmpereLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    void operator()(GridLayout& layout, auto&&... args) { core_type{layout}(args...); }

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



template<typename GridLayout, typename AMR_Types>
class OhmTransformer
{
    using info_type = core::OhmInfo;
    using core_type = core::Ohm<GridLayout>;
    using level_t   = AMR_Types::level_t;

public:
    explicit OhmTransformer(initializer::PHAREDict const& dict)
        : info_{info_type::FROM(dict)}
    {
    }

    void operator()(GridLayout& layout, auto&&... args) { core_type{info_, layout}(args...); }

    void operator()(level_t& level, auto& model, auto& B, auto& J, auto& E)
    {
        auto& electrons = model.state.electrons;
        auto& rm        = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, electrons, B, J, E))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            auto& n     = electrons.density();
            auto& Ve    = electrons.velocity();
            auto& Pe    = electrons.pressure();
            (*this)(layout, n, Ve, Pe, B, J, E);
        }
    }

    void operator()(level_t& level, auto& model, auto& B, auto& E)
    {
        (*this)(level, model, B, model.state.J, E);
    }

    info_type info_;
};


template<typename Model>
class OhmLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using info_type  = core::OhmInfo;
    using core_type  = core::Ohm<GridLayout>;

public:
    explicit OhmLevelTransformer(info_type const& info, level_t& level, Model& model)
        : info_{info}
        , level_{level}
        , model_{model}
    {
    }

    void operator()(GridLayout& layout, auto&&... args) { core_type{info_, layout}(args...); }

    void operator()(auto& B, auto& J, auto& E, auto& electrons)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, electrons, B, J, E))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            auto& n     = electrons.density();
            auto& Ve    = electrons.velocity();
            auto& Pe    = electrons.pressure();
            (*this)(layout, n, Ve, Pe, B, J, E);
        }
    }

    void operator()(auto& B, auto& E, auto& electrons) { (*this)(B, model_.state.J, E, electrons); }

    info_type info_;
    level_t& level_;
    Model& model_;
};

template<typename Model>
OhmLevelTransformer(core::OhmInfo, typename Model::amr_types::level_t&, Model&)
    -> OhmLevelTransformer<Model>;


} // namespace PHARE::solver



#endif /* PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP */
