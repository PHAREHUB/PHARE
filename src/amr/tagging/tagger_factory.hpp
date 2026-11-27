#ifndef PHARE_TAGGER_FACTORY_HPP
#define PHARE_TAGGER_FACTORY_HPP

#include "initializer/data_provider.hpp"

#include "core/def.hpp"
#include "core/utilities/meta/type_meta.hpp"

#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/tagging/concrete_tagger.hpp"
#include "amr/tagging/default_tagger_strategy.hpp"
#include "amr/tagging/gradient_scale_tagger_strategy.hpp"
#include "amr/tagging/lohner_cross_derivatives_tagger_strategy.hpp"
#include "amr/tagging/lohner_per_direction_tagger_strategy.hpp"
#include "amr/tagging/tagger.hpp"

#include "dict.hpp"

#include <memory>
#include <string>

namespace PHARE::amr
{
template<typename Model>
class TaggerFactory
{
public:
    TaggerFactory() = delete;

    NO_DISCARD static std::unique_ptr<Tagger<typename Model::amr_types>>
    make(PHARE::initializer::PHAREDict const& dict);
};

template<typename Model>
std::unique_ptr<Tagger<typename Model::amr_types>>
TaggerFactory<Model>::make(PHARE::initializer::PHAREDict const& dict)
{
    using HT  = ConcreteTagger<Model>;
    using HTS = DefaultTaggerStrategy<Model>;
    using LPD = LohnerPerDirectionTaggerStrategy<Model>;
    using LCD = LohnerCrossDerivativesTaggerStrategy<Model>;
    using GS  = GradientScaleTaggerStrategy<Model>;

    if constexpr (solver::is_hybrid_model_v<Model>)
    {
        std::string const methodName = dict["hybrid_method"];

        if (methodName == "default")
            return std::make_unique<HT>(std::make_unique<HTS>(dict));
        if (methodName == "lohner-per-direction")
            return std::make_unique<HT>(std::make_unique<LPD>(dict));
        if (methodName == "lohner-cross-derivatives")
            return std::make_unique<HT>(std::make_unique<LCD>(dict));
        if (methodName == "gradient-scale")
            return std::make_unique<HT>(std::make_unique<GS>(dict));
    }
    else if constexpr (solver::is_mhd_model_v<Model>)
    {
        std::string const methodName = dict["mhd_method"];

        if (methodName == "default")
            return std::make_unique<HT>(std::make_unique<HTS>(dict));
        if (methodName == "lohner-per-direction")
            return std::make_unique<HT>(std::make_unique<LPD>(dict));
        if (methodName == "lohner-cross-derivatives")
            return std::make_unique<HT>(std::make_unique<LCD>(dict));
        if (methodName == "gradient-scale")
            return std::make_unique<HT>(std::make_unique<GS>(dict));
    }
    else
        static_assert(core::dependent_false_v<Model>);

    return nullptr;
}


} // namespace PHARE::amr


#endif
