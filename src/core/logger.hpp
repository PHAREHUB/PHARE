#ifndef PHARE_CORE_LOGGER_HPP
#define PHARE_CORE_LOGGER_HPP

#include <cstdint>

#if !defined(PHARE_LOG_LEVEL)
#define PHARE_LOG_LEVEL 0 // 0 == off
#endif

namespace PHARE
{
constexpr static std::uint8_t LOG_LEVEL = PHARE_LOG_LEVEL;
}

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO)
#define PHARE_LOG_LINE_STR(str)                                                                    \
    std::cout << __FILE__ << ":" << __LINE__ << " - " << str << std::endl;
#else
#define PHARE_LOG_LINE_STR(str)
#endif
#define PHARE_LOG_LINE PHARE_LOG_LINE_STR("")

#if PHARE_WITH_CALIPER
#include "caliper/cali.h"

#define PHARE_LOG_START(lvl, str) CALI_MARK_BEGIN(str)
#define PHARE_LOG_STOP(lvl, str) CALI_MARK_END(str)
#define PHARE_LOG_SCOPE(lvl, str) PHARE::scope_log __phare_scope##__line__(lvl, str)

#else // !PHARE_WITH_CALIPER

#include "core/utilities/logger/logger_defaults.hpp"


#endif // PHARE_WITH_CALIPER

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

#endif /* PHARE_CORE_LOGGER_H */
