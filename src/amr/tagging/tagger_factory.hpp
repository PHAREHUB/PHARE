#ifndef PHARE_TAGGER_FACTORY_HPP
#define PHARE_TAGGER_FACTORY_HPP

#include <string>
#include <memory>

#include "tagger.hpp"
#include "hybrid_tagger.hpp"
#include "hybrid_tagger_strategy.hpp"
#include "default_hybrid_tagger_strategy.hpp"
#include "core/def.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::amr
{
template<typename PHARE_T>
class TaggerFactory
{
public:
    TaggerFactory() = delete;
    NO_DISCARD static std::unique_ptr<Tagger> make(PHARE::initializer::PHAREDict const& dict);
};

template<typename PHARE_T>
std::unique_ptr<Tagger> TaggerFactory<PHARE_T>::make(PHARE::initializer::PHAREDict const& dict)
{
    auto modelName  = dict["model"].template to<std::string>();
    auto methodName = dict["method"].template to<std::string>();

    if (modelName == "HybridModel")
    {
        using HybridModel = typename PHARE_T::HybridModel_t;
        using HT          = HybridTagger<HybridModel>;

        if (methodName == "default")
        {
            using HTS = DefaultHybridTaggerStrategy<HybridModel>;
            return std::make_unique<HT>(std::make_unique<HTS>(dict));
        }
    }
    return nullptr;
}


} // namespace PHARE::amr




#endif
