#ifndef PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP

#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "core/numerics/ohm/ohm.hpp"


namespace PHARE::solver
{



template<typename Model>
class FaradayLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;

    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i](); });
    }

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
        using field_type = VecField::field_type;

        if constexpr (core::is_field_tile_set_v<field_type>)
        {
            using Tile_vt = field_type::value_type::value_type;
            using V_t     = core::basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < B[0]().size(); ++tidx)
            {
                auto Bnw             = Bnew.template as<V_t>([&](auto& c) { return c()[tidx](); });
                auto const& tile_lay = B[0]()[tidx].layout();
                using TL             = std::remove_cvref_t<decltype(tile_lay)>;
                core::Faraday<TL>{tile_lay}(tt<V_t>(B, tidx), tt<V_t>(E, tidx), Bnw, dt);
            }
            for (std::uint8_t i = 0; i < 3; ++i)
                Bnew[i].sync_inner_ghosts();
        }
        else
        {
            core::Faraday<GridLayout>{layout}(B, E, Bnew, dt);
        }
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




template<typename Model>
class AmpereLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;

    template<typename V_t>
    V_t static tt(auto& vf, auto i)
    {
        return vf.template as<V_t>([&](auto& c) { return c()[i]; });
    }

public:
    explicit AmpereLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    template<typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField& J)
    {
        using field_type = VecField::field_type;

        if constexpr (core::is_field_tile_set_v<field_type>)
        {
            using Tile_vt = field_type::value_type;
            using V_t     = core::basic::TensorField<Tile_vt, 1>;

            for (std::size_t tidx = 0; tidx < J[0]().size(); ++tidx)
            {
                auto Jt              = J.template as<V_t>([&](auto& c) { return c()[tidx]; });
                auto const& tile_lay = J[0]()[tidx].layout();
                using TL             = std::remove_cvref_t<decltype(tile_lay)>;
                core::Ampere<TL>{tile_lay}(tt<V_t>(B, tidx), Jt);
            }
            for (std::uint8_t i = 0; i < 3; ++i)
                J[i].sync_inner_ghosts();
        }
        else
        {
            core::Ampere<GridLayout>{layout}(B, J);
        }

        core::check_tensor_field(J, layout);
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
