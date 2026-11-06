#ifndef PHARE_CORE_UTILITIES_ENV_HPP
#define PHARE_CORE_UTILITIES_ENV_HPP

#include "core/def.hpp"
#include "core/utilities/strings.hpp"

#include <string>
#include <cstdlib>
#include <optional>

namespace PHARE::core
{
NO_DISCARD inline std::optional<std::string> get_env(std::string const& key)
{
    if (char const* val = std::getenv(key.c_str()))
        return std::string{val};
    return std::nullopt;
}

NO_DISCARD inline std::string get_env(std::string const& key, std::string const& _default)
{
    if (auto e = get_env(key))
        return *e;
    return _default;
}


template<typename T>
NO_DISCARD inline T get_env_as(std::string const& key, T const& t)
{
    if (auto e = get_env(key))
        return from_string<T>(*e);
    return t;
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_ENV_HPP
