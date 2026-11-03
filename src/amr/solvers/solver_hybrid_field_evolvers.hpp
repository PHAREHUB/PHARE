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
        amr::onLevel(amr::make_model_level_accessor(level_, model_, electrons, B, J, E),
                     [&](auto& view) {
                         auto& [electrons_v, B_v, J_v, E_v] = view.args;
                         (*this)(view.layout, electrons_v.density(), electrons_v.velocity(),
                                 electrons_v.pressure(), B_v, J_v, E_v);
                     });
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
