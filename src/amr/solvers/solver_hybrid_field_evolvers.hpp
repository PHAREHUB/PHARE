#ifndef PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/data/field/field_scratch.hpp"
#include "solver_field_evolvers.hpp"

namespace PHARE::solver
{


template<typename Model>
class OhmLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using VecField   = Model::vecfield_type;
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
        auto& rm          = *model_.resourcesManager;
        auto pressureTerm = VecField{"pressureTerm", core::HybridQuantity::Vector::E};
        for (auto& patch : rm.enumerate(level_, electrons, B, J, E))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            // the electron closure writes the pressure term of this patch in the shared field
            // scratch, Ohm reads it; the scratch is held meanwhile, so neither may take it
            core::ScratchWrite scratch;
            setOnScratch_(pressureTerm, layout, scratch.get());
            electrons.computePressureTerm(layout, pressureTerm);
            (*this)(layout, electrons.density(), electrons.velocity(), pressureTerm, B, J, E);
            pressureTerm.setBuffer(nullptr);
        }
    }

    void operator()(auto& B, auto& E, auto& electrons) { (*this)(B, model_.state.J, E, electrons); }

    info_type info_;
    level_t& level_;
    Model& model_;

private:
    // each component of v gets its share of the scratch, shaped for this layout
    static void setOnScratch_(VecField& v, GridLayout const& layout, auto& buffer)
    {
        auto const shapes = core::for_N<VecField::size(), core::for_N_R_mode::make_array>(
            [&](auto i) { return layout.allocSize(v[i].physicalQuantity()); });
        auto total = std::size_t{0};
        for (auto const& shape : shapes)
            total += core::product(shape, std::size_t{1});

        auto* data = buffer.get_no_copy(total).data();
        for (std::size_t i = 0; i < shapes.size(); ++i)
        {
            v[i] = typename VecField::field_type{v[i].name(), v[i].physicalQuantity(), data,
                                                 shapes[i]};
            data += core::product(shapes[i], std::size_t{1});
        }
    }
};

template<typename Model>
OhmLevelTransformer(core::OhmInfo, typename Model::amr_types::level_t&, Model&)
    -> OhmLevelTransformer<Model>;


} // namespace PHARE::solver


#endif /* PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP */
