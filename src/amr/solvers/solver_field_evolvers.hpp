#ifndef PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "amr/resources_manager/amr_utils.hpp"


namespace PHARE::solver
{


template<typename GridLayout, typename AMR_Types>
class FaradayTransformer
{
    using core_type = core::Faraday<GridLayout>;
    using level_t   = AMR_Types::level_t;

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

template<typename GridLayout, typename AMR_Types>
class AmpereTransformer
{
    using core_type = core::Ampere<GridLayout>;
    using level_t   = AMR_Types::level_t;

public:
    void operator()(GridLayout& layout, auto&&... args) { core_type{layout}(args...); }

    void operator()(level_t& level, auto& model, auto& B, auto& J)
    {
        auto& rm    = *model.resourcesManager;
        auto& state = model.state;
        for (auto& patch : rm.enumerate(level, state, B, J))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            (*this)(layout, B, J);
        }
    }
};


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



} // namespace PHARE::solver



#endif /* PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP */
