#ifndef PHARE_LEVEL_INITIALIZER_FACTORY_H
#define PHARE_LEVEL_INITIALIZER_FACTORY_H

#include "hybrid_level_initializer.hpp"
#include "level_initializer.hpp"
#include "initializer/data_provider.hpp"

#include <memory>
#include <string>

namespace PHARE
{
namespace solver
{
    template<typename HybridModel>
    class LevelInitializerFactory
    {
        using AMRTypes = typename HybridModel::amr_types;

    public:
        static std::unique_ptr<LevelInitializer<AMRTypes>>
        create(std::string modelName, PHARE::initializer::PHAREDict const& dict)
        {
            if (modelName == "HybridModel")
            {
                return std::make_unique<HybridLevelInitializer<HybridModel>>(dict);
            }
            return nullptr;
        }
    };

} // namespace solver
} // namespace PHARE



#endif
