#ifndef PHARE_CORE_LOGGER_HPP
#define PHARE_CORE_LOGGER_HPP

#include <cstdint>
#include <unordered_map>

#if !defined(PHARE_LOG_LEVEL)
#define PHARE_LOG_LEVEL 1 // 0 == off
#endif

namespace PHARE
{

constexpr static std::uint8_t LOG_LEVEL = PHARE_LOG_LEVEL;

}

// clang-format off
#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO) || defined(PHARE_FORCE_LOG_LINE)
#include <sstream>
#include <iostream>
#define PHARE_LOG_LINE_STR(str)                                                                    \
    std::cout << __FILE__ << ":" << __LINE__ << " - " << str << std::endl;
#define PHARE_LOG_LINE_SS(s) PHARE_LOG_LINE_STR((std::stringstream{} << s).str())

#else // LOGGING DISABLED

#define PHARE_LOG_LINE_STR(str) {} // allow if/for/etc with no scope
#define PHARE_LOG_LINE_SS(str)  {} // allow if/for/etc with no scope

#endif // LOGGING


#if PHARE_WITH_CALIPER
#include "caliper/cali.h"

#define PHARE_LOG_START(lvl, str) CALI_MARK_BEGIN(str)
#define PHARE_LOG_STOP(lvl, str) CALI_MARK_END(str)
#define PHARE_LOG_SCOPE(lvl, str) PHARE::scope_log __phare_scope##__line__(lvl, str)

#else // !PHARE_WITH_CALIPER

#include "core/utilities/logger/logger_defaults.hpp"
#endif // PHARE_WITH_CALIPER

// clang-format on

#include <string>
#include <utility>

namespace PHARE
{
struct scope_log
{
    scope_log(int&& i_, std::string&& str)
        : i{i_}
        , key{std::move(str)}
    {
        if (i <= LOG_LEVEL)
        {
            PHARE_LOG_START(i, key.c_str());
        }
    }
    ~scope_log()
    {
        if (i <= LOG_LEVEL)
        {
            PHARE_LOG_STOP(i, key.c_str());
        }
    }

    int i;
    std::string key;
};
} // namespace PHARE


// FUNCTION COUNT LOGGING
namespace PHARE::core
{

namespace detail
{
#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO) || defined(PHARE_FORCE_FUNC_COUNT)
    bool static constexpr inline log_counting_active = 1;
#else
    bool static constexpr inline log_counting_active = 0;
#endif
} // namespace detail

struct FunctionCountMonitor
{
    ~FunctionCountMonitor()
    {
        if (detail::log_counting_active)
            for ([[maybe_unused]] auto const& [k, v] : ops)
            {
                PHARE_LOG_LINE_SS(k << " " << v);
            }
    }

    void operator()(std::string const& l) { ops[s + "::" + l] += 1; }

    std::string s;
    std::unordered_map<std::string, std::size_t> ops{};
};

} // namespace PHARE::core

// #if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO) || defined(PHARE_FORCE_FUNC_COUNT)
// #define PHARE_FUNC_COUNT(instance, string) instance(string);
// #else
#define PHARE_FUNC_COUNT(a, b) // noop
// #endif


#endif /* PHARE_CORE_LOGGER_H */
