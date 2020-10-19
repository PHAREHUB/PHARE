
#ifndef PHARE_TAGGER_FACTORY_H
#define PHARE_TAGGER_FACTORY_H

#include <string>
#include <memory>

#include "tagger.h"
#include "hybrid_tagger.h"
#include "hybrid_tagger_strategy.h"
#include "scaledavg_hybrid_tagger_strategy.h"

namespace PHARE::amr
{
template<typename PHARE_T>
class TaggerFactory
{
public:
    TaggerFactory() = delete;
    static std::unique_ptr<Tagger> make(std::string modelName, std::string methodName);
};

template<typename PHARE_T>
std::unique_ptr<Tagger> TaggerFactory<PHARE_T>::make(std::string modelName, std::string methodName)
{
    if (modelName == "HybridModel")
    {
        using HybridModel = typename PHARE_T::HybridModel_t;
        using HT          = HybridTagger<HybridModel>;
        using HTS         = ScaledAvgHybridTaggerStrategy<HybridModel>;
        return std::make_unique<HT>(std::make_unique<HTS>());
    }
    return nullptr;
}


} // namespace PHARE::amr




#endif
