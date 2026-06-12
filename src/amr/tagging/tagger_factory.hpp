#ifndef PHARE_TAGGER_FACTORY_HPP
#define PHARE_TAGGER_FACTORY_HPP

#include "core/def.hpp"

#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"

#include "initializer/data_provider.hpp"

#include "tagger.hpp"
#include "concrete_tagger.hpp"
#include "default_tagger_strategy.hpp"

#include <string>
#include <memory>

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
    using HT  = ConcreteTagger<Model>;
    using HTS = DefaultTaggerStrategy<Model>;

    if constexpr (solver::is_hybrid_model_v<Model>)
    {
        std::string const methodName = dict["hybrid_method"];

        if (methodName == "default")
            return std::make_unique<HT>(std::make_unique<HTS>(dict));
    }
    else if constexpr (solver::is_mhd_model_v<Model>)
    {
        std::string const methodName = dict["mhd_method"];

        if (methodName == "default")
            return std::make_unique<HT>(std::make_unique<HTS>(dict));
    }
    else
        static_assert(core::dependent_false_v<Model>);

    return nullptr;
}


} // namespace PHARE::amr


#endif
