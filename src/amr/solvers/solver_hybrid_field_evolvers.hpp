#ifndef PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP

#include "core/numerics/ohm/ohm.hpp"
#include "core/data/field/field_tiles.hpp"

#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/resources_manager/amr_utils.hpp"

namespace PHARE::solver
{

class OhmSingleTransformer
{
    using info_type = core::OhmInfo;

    template<typename GridLayout>
    void operate(GridLayout const& layout, auto&&... args)
    {
        core::Ohm<GridLayout>{info_, layout}(args...);
    }

public:
    explicit OhmSingleTransformer(info_type const& info)
        : info_{info}
    {
    }

    template<typename GridLayout, typename VecField, typename Field>
    void operator()(GridLayout const& layout, Field const& n, VecField const& Ve, Field const& Pe,
                    VecField const& B, VecField const& J, VecField& Enew)
        requires(core::is_field_v<Field>)
    {
        operate(layout, n, Ve, Pe, B, J, Enew);
    }

    template<typename GridLayout, typename VecField, typename Field>
    void operator()(GridLayout const& /*layout*/, Field const& n, VecField const& Ve,
                    Field const& Pe, VecField const& B, VecField const& J, VecField& Enew)
        requires(core::is_field_tile_set_v<Field>)
    {
        for (std::size_t i = 0; i < core::tile_count(Enew); ++i)
            on_tile(i, n, Ve, Pe, B, J, Enew);

        core::sync_inner_ghosts(Enew);
    }

    // one tile's worth of work, Enew ghosts are not synced
    void on_tile(std::size_t const tile_idx, auto const& n, auto const& Ve, auto const& Pe,
                 auto const& B, auto const& J, auto& Enew)
    {
        core::tile_exec_with_layout_at(
            tile_idx, [&](auto& layout, auto&&... args) { operate(layout, args...); }, n, Ve, Pe, B,
            J, Enew);
    }

    info_type info_;
};

template<typename Model>
class OhmLevelTransformer : public OhmSingleTransformer
{
    using Super      = OhmSingleTransformer;
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using info_type  = core::OhmInfo;

public:
    explicit OhmLevelTransformer(info_type const& info, level_t& level, Model& model)
        : Super{info}
        , level_{level}
        , model_{model}
    {
    }

    template<typename VecField>
    void operator()(VecField& B, VecField& J, VecField& E, auto& electrons)
        requires(not core::has_tiled_field_type_c<VecField>)
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

    template<typename VecField>
    void operator()(VecField& B, VecField& J, VecField& E, auto& electrons)
        requires(core::has_tiled_field_type_c<VecField>)
    {
        auto accessor = amr::make_model_level_accessor(level_, model_, electrons, B, J, E);
        tiled_level_transform(
            accessor, [](auto& view) -> auto& { return std::get<3>(view.args); },
            [this](auto& view, auto const tile_idx) {
                auto& [electrons_v, B_v, J_v, E_v] = view.args;
                Super::on_tile(tile_idx, electrons_v.density(), electrons_v.velocity(),
                               electrons_v.pressure(), B_v, J_v, E_v);
            });
    }

    void operator()(auto& B, auto& E, auto& electrons) { (*this)(B, model_.state.J, E, electrons); }

    level_t& level_;
    Model& model_;
};

template<typename Model>
OhmLevelTransformer(core::OhmInfo, typename Model::amr_types::level_t&, Model&)
    -> OhmLevelTransformer<Model>;


// Ve and Pe from the ions moments and J, per patch or per tile on the thread pools
template<typename Model>
class ElectronsLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using field_type = Model::field_type;

public:
    explicit ElectronsLevelTransformer(level_t& level, Model& model)
        : level_{level}
        , model_{model}
    {
    }

    void operator()(auto& electrons)
        requires(not core::is_field_tile_set_v<field_type>)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, electrons))
            electrons.update(amr::layoutFromPatch<GridLayout>(*patch));
    }

    void operator()(auto& electrons)
        requires(core::is_field_tile_set_v<field_type>)
    {
        auto accessor = amr::make_model_level_accessor(level_, model_, electrons);
        tiled_level_exec(
            accessor, [](auto& view) -> auto& { return std::get<0>(view.args).density(); },
            [](auto& view, auto const tile_idx) { std::get<0>(view.args).update(tile_idx); });
    }

    level_t& level_;
    Model& model_;
};

template<typename Model>
ElectronsLevelTransformer(typename Model::amr_types::level_t&, Model&)
    -> ElectronsLevelTransformer<Model>;

} // namespace PHARE::solver

#endif /* PHARE_AMR_SOLVERS_SOLVER_HYBRID_FIELD_EVOLVERS_HPP */
