
#ifndef PHARE_HYBRID_TAGGER_HPP
#define PHARE_HYBRID_TAGGER_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "tagger.hpp"
#include "hybrid_tagger_strategy.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/types/amr_types.hpp"

#include <SAMRAI/pdat/CellData.h>

#include <memory>
#include <utility>
#include <stdexcept>




namespace PHARE::amr
{
template<typename HybridModel>
class HybridTagger : public Tagger
{
    using patch_t         = typename Tagger::patch_t;
    using amr_t           = PHARE::amr::SAMRAI_Types;
    using IPhysicalModel  = PHARE::solver::IPhysicalModel<amr_t>;
    using gridlayout_type = typename HybridModel::gridlayout_type;


public:
    HybridTagger(std::unique_ptr<HybridTaggerStrategy<HybridModel>> strat)
        : Tagger{"HybridTagger"}
        , strat_{std::move(strat)}
    {
    }

    void tag(IPhysicalModel& model, patch_t& patch, int tag_index) override;

private:
    std::unique_ptr<HybridTaggerStrategy<HybridModel>> strat_;
};




//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------




template<typename HybridModel>
void HybridTagger<HybridModel>::tag(PHARE::solver::IPhysicalModel<amr_t>& model, patch_t& patch,
                                    int tag_index)
{
    if (strat_)
    {
        auto& hybridModel   = dynamic_cast<HybridModel&>(model);
        auto layout         = PHARE::amr::layoutFromPatch<gridlayout_type>(patch);
        auto modelIsOnPatch = hybridModel.setOnPatch(patch);
        auto pd   = dynamic_cast<SAMRAI::pdat::CellData<int>*>(patch.getPatchData(tag_index).get());
        auto tags = pd->getPointer();
        strat_->tag(hybridModel, layout, tags);


        // These tags will be saved even if they are not used in diags during this advance
        // hybridModel.tags may contain vectors for patches and levels that no longer exist
        auto key = std::to_string(patch.getPatchLevelNumber()) + "_"
                   + amr::to_string(patch.getGlobalId());

        auto nCells = core::product(layout.nbrCells());

        bool item_exists_and_valid
            = hybridModel.tags.count(key) and hybridModel.tags[key]->size() == nCells;

        if (!item_exists_and_valid)
        {
            using Map_value_type = typename std::decay_t<decltype(hybridModel.tags)>::mapped_type;


            hybridModel.tags[key]
                = std::make_shared<typename Map_value_type::element_type>(layout.nbrCells());
        }

        auto nbrCells = layout.nbrCells();
        auto tagsv = core::NdArrayView<HybridModel::dimension, int>(hybridModel.tags[key]->data(),
                                                                    layout.nbrCells());
        auto tagsvF
            = core::NdArrayView<HybridModel::dimension, int, false>(tags, layout.nbrCells());
        if constexpr (HybridModel::dimension == 2)
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
