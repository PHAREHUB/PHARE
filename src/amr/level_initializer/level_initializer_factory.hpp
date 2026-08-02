#ifndef PHARE_LEVEL_INITIALIZER_FACTORY_HPP
#define PHARE_LEVEL_INITIALIZER_FACTORY_HPP

#include "mhd_level_initializer.hpp"
#include "hybrid_level_initializer.hpp"
#include "level_initializer.hpp"
#include "initializer/data_provider.hpp"
#include "core/def.hpp"

#include <memory>
#include <string>
#include <type_traits>

namespace PHARE
{
namespace solver
{
    // Variadic over the enabled initializers only — an MHD-only build never names
    // HybridLevelInitializer<HybridModel>, so no hybrid type is instantiated for it.
    template<typename AMRTypes, typename... Initializers>
    class LevelInitializerFactory
    {
    public:
        NO_DISCARD static std::unique_ptr<LevelInitializer<AMRTypes>>
        create(std::string modelName, PHARE::initializer::PHAREDict const& dict)
        {
            std::unique_ptr<LevelInitializer<AMRTypes>> result;
            ((result = tryCreate<Initializers>(modelName, dict)) || ...);
            return result;
        }

    private:
        template<typename Initializer>
        static std::unique_ptr<LevelInitializer<AMRTypes>>
        tryCreate(std::string const& modelName, PHARE::initializer::PHAREDict const& dict)
        {
            if (modelName != Initializer::model_type::model_name)
                return {};
            if constexpr (std::is_constructible_v<Initializer,
                                                   PHARE::initializer::PHAREDict const&>)
                return std::make_unique<Initializer>(dict);
            else
                return std::make_unique<Initializer>();
        }
    };

} // namespace solver
} // namespace PHARE



#endif
