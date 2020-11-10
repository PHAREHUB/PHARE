
#ifndef PHARE_HYBRID_TAGGER_H
#define PHARE_HYBRID_TAGGER_H

#include "tagger.h"
#include "solver/physical_models/hybrid_model.h"
#include "amr/types/amr_types.h"

#include <SAMRAI/pdat/CellData.h>

#include <memory>
#include <utility>
#include <stdexcept>




namespace PHARE::amr
{
template<typename HybridModel>
class HybridTaggerStrategy
{
    using gridlayout_type = typename HybridModel::gridlayout_type;

public:
    virtual void tag(HybridModel& model, gridlayout_type const& layout, int* tags) const = 0;
    virtual ~HybridTaggerStrategy()                                                      = 0;
};

template<typename HybridModel>
HybridTaggerStrategy<HybridModel>::~HybridTaggerStrategy()
{
}


template<typename HybridModel>
class ScaledAvgHybridTaggerStrategy : public HybridTaggerStrategy<HybridModel>
{
    using gridlayout_type = typename HybridModel::gridlayout_type;

public:
    void tag(HybridModel& model, gridlayout_type const& layout, double* tags) const override;
};

template<typename HybridModel>
void ScaledAvgHybridTaggerStrategy<HybridModel>::tag(HybridModel& model,
                                                     gridlayout_type const& layout,
                                                     double* tags) const
{
}




template<typename HybridModel>
class HybridTagger : public Tagger
{
    using patch_t         = typename Tagger::patch_t;
    using amr_t           = PHARE::amr::SAMRAI_Types;
    using gridlayout_type = typename HybridModel::gridlayout_type;


public:
    HybridTagger(std::unique_ptr<HybridTaggerStrategy<HybridModel>> strat)
        : strat_{std::move(strat)}
    {
    }

    void tag(PHARE::solver::IPhysicalModel<amr_t>& model, patch_t& patch, int tag_index) override;

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
    }
    else
        throw std::runtime_error("invalid tagging strategy");
}

} // namespace PHARE::amr

#endif
