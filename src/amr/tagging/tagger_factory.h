
#ifndef PHARE_TAGGER_FACTORY_H
#define PHARE_TAGGER_FACTORY_H

#include <string>
#include <memory>

#include "tagger.h"
#include "hybrid_tagger.h"

namespace PHARE::amr
{
class TaggerFactory
{
public:
    static std::unique_ptr<Tagger> make(std::string modelName, std::string methodName);
};


std::unique_ptr<Tagger> TaggerFactory::make(std::string modelName, std::string methodName)
{
    if (modelName == "hybridModel")
    {
        return /*std::make_unique<HybridTagger>(std::make_unique<ScaledAvgHybridTaggerStrategy>())*/
            nullptr;
    }
    return nullptr;
}


} // namespace PHARE::amr




#endif
