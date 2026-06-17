#ifndef PHARE_CONCRETE_TAGGER_HPP
#define PHARE_CONCRETE_TAGGER_HPP


#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/utilities/types.hpp"

#include "tagger.hpp"
#include "tagging_criteria.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/types/amr_types.hpp"

#include "initializer/data_provider.hpp"

#include <SAMRAI/pdat/CellData.h>

#include <memory>
#include <string>
#include <vector>
#include <stdexcept>




namespace PHARE::amr
{
template<typename Model>
class ConcreteTagger : public Tagger
{
    using patch_t         = typename Tagger::patch_t;
    using amr_t           = PHARE::amr::SAMRAI_Types;
    using IPhysicalModel  = PHARE::solver::IPhysicalModel<amr_t>;
    using gridlayout_type = typename Model::gridlayout_type;

    static auto constexpr dimension = Model::dimension;

    // hybrid stores B in electromag, MHD stores it directly in the state; static
    // so it can also seed the field_type alias below.
    static auto& modelB(Model& model)
    {
        if constexpr (solver::is_hybrid_model_v<Model>)
            return model.state.electromag.B;
        else if constexpr (solver::is_mhd_model_v<Model>)
            return model.state.B;
        else
            static_assert(core::dependent_false_v<Model>);
    }

    using vecfield_t = std::decay_t<decltype(modelB(std::declval<Model&>()))>;
    using field_type = typename vecfield_t::field_type;

    struct QuantityTag
    {
        std::string name;
        double threshold;
    };

public:
    ConcreteTagger(initializer::PHAREDict const& dict)
        : Tagger{Model::model_name == "HybridModel" ? "HybridTagger" : "MHDTagger"}
        , method_{parseTaggingMethod(dict["method"].template to<std::string>())}
    {
        auto const nbrQuantities = dict["nbr_quantities"].template to<int>();
        for (int i = 0; i < nbrQuantities; ++i)
        {
            auto const path        = "Q" + std::to_string(i);
            std::string const name = dict[path]["name"].template to<std::string>();
            double const threshold = dict[path]["threshold"].template to<double>();
            if (!supported_(name))
                throw std::runtime_error("tagging quantity '" + name
                                         + "' is not supported for this model");
            quantities_.push_back({name, threshold});
        }
    }

    void tag(IPhysicalModel& model, patch_t& patch, int tag_index) override;

private:
    void tagCells_(Model& model, gridlayout_type const& layout, int* tags) const;
    std::vector<field_type const*> components_(Model& model, std::string const& qty) const;
    static bool supported_(std::string const& qty);

    TaggingMethod method_;
    std::vector<QuantityTag> quantities_;
};




//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------




template<typename Model>
bool ConcreteTagger<Model>::supported_(std::string const& qty)
{
    if (qty == "B")
        return true;
    if (qty == "rho")
        return solver::is_mhd_model_v<Model>;
    return false;
}


template<typename Model>
std::vector<typename ConcreteTagger<Model>::field_type const*>
ConcreteTagger<Model>::components_(Model& model, std::string const& qty) const
{
    if (qty == "B")
    {
        auto&& [bx, by, bz] = modelB(model)();
        return {&bx, &by, &bz};
    }
    if (qty == "rho")
    {
        if constexpr (solver::is_mhd_model_v<Model>)
            return {&model.state.rho};
        else
            throw std::runtime_error("rho tagging not yet supported on hybrid");
    }
    throw std::runtime_error("unknown tagging quantity '" + qty + "'");
}


template<typename Model>
void ConcreteTagger<Model>::tagCells_(Model& model, gridlayout_type const& layout, int* tags) const
{
    auto const nbrCells = layout.nbrCells();

    // SAMRAI tags int* buffer is FORTRAN ordering so we set false to the view; it
    // has no ghost cells (one entry per physical cell).
    bool constexpr c_ordering = false;
    auto tagsv                = core::NdArrayView<dimension, int, c_ordering>(tags, nbrCells);

    // union over quantities: a cell is tagged if ANY quantity's indicator exceeds
    // its own threshold. Zero once, then only ever set to 1.
    std::fill(tags, tags + core::product(nbrCells), 0);

    // We loop on cell indexes for all quantities regardless of their centering
    // (co-indexing, see tagging_criteria.hpp). The centered stencil reaches +/-1,
    // valid at every physical cell since nbrGhosts >= 1, so no clamp is needed.
    std::array<std::uint32_t, dimension> start;
    start[0] = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::X);
    if constexpr (dimension > 1)
        start[1] = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::Y);
    if constexpr (dimension > 2)
        start[2] = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::Z);

    for (auto const& q : quantities_)
    {
        auto const comps = components_(model, q.name);

        auto const indicator = [&](std::array<std::uint32_t, dimension> const& idx) {
            return method_ == TaggingMethod::Lohner ? lohnerIndicator<dimension>(comps, idx)
                                                    : defaultIndicator<dimension>(comps, idx);
        };

        if constexpr (dimension == 1)
        {
            for (std::uint32_t iCell = 0; iCell < nbrCells[0]; ++iCell)
                if (indicator({start[0] + iCell}) > q.threshold)
                    tagsv(iCell) = 1;
        }
        else if constexpr (dimension == 2)
        {
            for (std::uint32_t ix = 0; ix < nbrCells[0]; ++ix)
                for (std::uint32_t iy = 0; iy < nbrCells[1]; ++iy)
                    if (indicator({start[0] + ix, start[1] + iy}) > q.threshold)
                        tagsv(ix, iy) = 1;
        }
        else if constexpr (dimension == 3)
        {
            for (std::uint32_t ix = 0; ix < nbrCells[0]; ++ix)
                for (std::uint32_t iy = 0; iy < nbrCells[1]; ++iy)
                    for (std::uint32_t iz = 0; iz < nbrCells[2]; ++iz)
                        if (indicator({start[0] + ix, start[1] + iy, start[2] + iz}) > q.threshold)
                            tagsv(ix, iy, iz) = 1;
        }
    }
}


template<typename Model>
void ConcreteTagger<Model>::tag(PHARE::solver::IPhysicalModel<amr_t>& model, patch_t& patch,
                                int tag_index)
{
    auto& concreteModel = dynamic_cast<Model&>(model);
    auto layout         = PHARE::amr::layoutFromPatch<gridlayout_type>(patch);
    auto modelIsOnPatch = concreteModel.setOnPatch(patch);
    auto pd   = dynamic_cast<SAMRAI::pdat::CellData<int>*>(patch.getPatchData(tag_index).get());
    auto tags = pd->getPointer();
    tagCells_(concreteModel, layout, tags);


    // These tags will be saved even if they are not used in diags during this advance
    // concreteModel.tags may contain vectors for patches and levels that no longer exist
    auto key = std::to_string(patch.getPatchLevelNumber()) + "_"
               + amr::to_string(patch.getGlobalId());

    auto nCells = core::product(layout.nbrCells());

    bool item_exists_and_valid
        = concreteModel.tags.count(key) and concreteModel.tags[key]->size() == nCells;

    if (!item_exists_and_valid)
    {
        using Map_value_type = typename std::decay_t<decltype(concreteModel.tags)>::mapped_type;


        concreteModel.tags[key]
            = std::make_shared<typename Map_value_type::element_type>(layout.nbrCells());
    }

    auto nbrCells = layout.nbrCells();
    auto tagsv    = core::NdArrayView<Model::dimension, int>(concreteModel.tags[key]->data(),
                                                             layout.nbrCells());
    auto tagsvF   = core::NdArrayView<Model::dimension, int, false>(tags, layout.nbrCells());
    if constexpr (Model::dimension == 2)
    {
        for (auto iTag_x = 0u; iTag_x < nbrCells[0]; ++iTag_x)
        {
            for (auto iTag_y = 0u; iTag_y < nbrCells[1]; ++iTag_y)
            {
                tagsv(iTag_x, iTag_y) = tagsvF(iTag_x, iTag_y);
            }
        }
    }
}

} // namespace PHARE::amr

#endif
