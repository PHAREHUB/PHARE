#ifndef PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "solver_field_evolvers.hpp"

namespace PHARE::solver
{


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
        auto& rm   = *model_.resourcesManager;
        auto& ions = model_.state.ions;
        for (auto& patch : rm.enumerate(level_, electrons, ions, B, J, E))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            // Ohm takes the electron flux n Ve = Ji - J on the E nodes, formed where J and E
            // live (#1577), not Ve on the moment nodes
            (*this)(layout, ions.chargeDensity(), electrons.flux(), electrons.pressureTerm(), B, J,
                    E);
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


#endif /* PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP */
