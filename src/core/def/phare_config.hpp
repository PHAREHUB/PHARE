#ifndef PHARE_CORE_DEF_PHARE_CONFIG_HPP
#define PHARE_CORE_DEF_PHARE_CONFIG_HPP

#include <string>


#if __has_include("core/def/_gen_sys.hpp")
#include "core/def/_gen_sys.hpp"
#else

namespace PHARE
{
std::unordered_map<std::string, std::string> build_config()
{
    return {{"PHARE_CONFIG_ERROR", "PHARE NOT BUILT WITH CONFIGURATOR!"}};
}
} // namespace PHARE

#endif

#endif /*PHARE_CORE_DEF_PHARE_CONFIG_HPP*/
