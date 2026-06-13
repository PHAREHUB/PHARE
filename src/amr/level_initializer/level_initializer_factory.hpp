#ifndef PHARE_LEVEL_INITIALIZER_FACTORY_HPP
#define PHARE_LEVEL_INITIALIZER_FACTORY_HPP

#include "mhd_level_initializer.hpp"
#include "hybrid_level_initializer.hpp"
#include "level_initializer.hpp"
#include "initializer/data_provider.hpp"
#include "core/def.hpp"

#include <memory>
#include <string>

namespace PHARE
{
namespace solver
{
    template<typename HybridModel, typename MHDModel>
    class LevelInitializerFactory
    {
        using AMRTypes = typename HybridModel::amr_types;

    public:
        NO_DISCARD static std::unique_ptr<LevelInitializer<AMRTypes>>
        create(std::string modelName, PHARE::initializer::PHAREDict const& dict)
        {
            if (modelName == "HybridModel")
            {
                return std::make_unique<HybridLevelInitializer<HybridModel>>(dict);
            }
            else if (modelName == "MHDModel")
            {
                return std::make_unique<MHDLevelInitializer<MHDModel>>();
            }
            return nullptr;
        }
    };

} // namespace solver
} // namespace PHARE



#endif
