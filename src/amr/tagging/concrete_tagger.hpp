
#ifndef PHARE_HYBRID_TAGGER_HPP
#define PHARE_HYBRID_TAGGER_HPP


#include "core/def/phare_mpi.hpp"

#include "tagger.hpp"
#include "tagger_strategy.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/types/amr_types.hpp"

#include <SAMRAI/pdat/CellData.h>

#include <memory>
#include <utility>
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

public:
    ConcreteTagger(std::unique_ptr<TaggerStrategy<Model>> strat)
        : Tagger{Model::model_name == "HybridModel" ? "HybridTagger" : "MHDTagger"}
        , strat_{std::move(strat)}
    {
    }

    void tag(IPhysicalModel& model, patch_t& patch, int tag_index);

private:
    std::unique_ptr<TaggerStrategy<Model>> strat_;
};




//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------




template<typename Model>
void ConcreteTagger<Model>::tag(PHARE::solver::IPhysicalModel<amr_t>& model, patch_t& patch,
                                int tag_index)
{
    if (strat_)
    {
        auto& concreteModel = dynamic_cast<Model&>(model);
        auto layout         = PHARE::amr::layoutFromPatch<gridlayout_type>(patch);
        auto modelIsOnPatch = concreteModel.setOnPatch(patch);
        auto pd   = dynamic_cast<SAMRAI::pdat::CellData<int>*>(patch.getPatchData(tag_index).get());
        auto tags = pd->getPointer();
        strat_->tag(concreteModel, layout, tags);


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
    else
        throw std::runtime_error("invalid tagging strategy");
}

} // namespace PHARE::amr

#endif
