#ifndef PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "solver_field_evolvers.hpp"

namespace PHARE::solver
{


template<typename Model>
class OhmLevelTransformer : public core::OhmSingleTransformer
{
    using Super      = OhmSingleTransformer;
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using info_type  = core::OhmInfo;

    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i](); });
    }

public:
    explicit OhmLevelTransformer(info_type const& info, level_t& level, Model& model)
        : Super{info}
        , level_{level}
        , model_{model}
    {
    }

    void operator()(auto& B, auto& J, auto& E, auto& electrons)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, electrons, B, J, E))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            auto& n     = electrons.density();
            auto& Ve    = electrons.velocity();
            auto& Pe    = electrons.pressure();
            Super::operator()(layout, n, Ve, Pe, B, J, E);
        }
    }

    void operator()(auto& B, auto& E, auto& electrons) { (*this)(B, model_.state.J, E, electrons); }

    level_t& level_;
    Model& model_;
};

template<typename Model>
OhmLevelTransformer(core::OhmInfo, typename Model::amr_types::level_t&, Model&)
    -> OhmLevelTransformer<Model>;

} // namespace PHARE::solver


#endif /* PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP */
