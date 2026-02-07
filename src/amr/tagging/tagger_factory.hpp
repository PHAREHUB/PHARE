#ifndef PHARE_TAGGER_FACTORY_HPP
#define PHARE_TAGGER_FACTORY_HPP

#include <string>
#include <memory>

#include "tagger.hpp"
#include "concrete_tagger.hpp"
#include "tagger_strategy.hpp"
#include "default_tagger_strategy.hpp"
#include "core/def.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::amr
{
template<typename Model>
class TaggerFactory
{
public:
    TaggerFactory() = delete;
    NO_DISCARD static std::unique_ptr<Tagger> make(PHARE::initializer::PHAREDict const& dict);
};

template<typename Model>
std::unique_ptr<Tagger> TaggerFactory<Model>::make(PHARE::initializer::PHAREDict const& dict)
{
    auto modelName = Model::model_name;

    if (modelName == "HybridModel")
    {
        auto methodName = dict["hybrid_method"].template to<std::string>();
        using HT        = ConcreteTagger<Model>;

        if (methodName == "default")
        {
            using HTS = DefaultTaggerStrategy<Model>;
            return std::make_unique<HT>(std::make_unique<HTS>(dict));
        }
    }
    else if (modelName == "MHDModel")
    {
        auto methodName = dict["mhd_method"].template to<std::string>();
        using HT        = ConcreteTagger<Model>;

        if (methodName == "default")
        {
            using HTS = DefaultTaggerStrategy<Model>;
            return std::make_unique<HT>(std::make_unique<HTS>(dict));
        }
    }
    return nullptr;
}


} // namespace PHARE::amr




#endif
