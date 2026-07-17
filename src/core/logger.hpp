#ifndef PHARE_CORE_LOGGER_HPP
#define PHARE_CORE_LOGGER_HPP

#include "core/def.hpp" // IWYU pragma: keep

#include <chrono>
#include <string>
#include <cstdint>
#include <utility>

#if !defined(PHARE_LOG_LEVEL)
#define PHARE_LOG_LEVEL 2 // 0 == off
#endif

namespace PHARE
{
constexpr static std::uint8_t LOG_LEVEL = PHARE_LOG_LEVEL;
} // namespace PHARE

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO) || defined(PHARE_FORCE_LOG_LINE)
#include <sstream>  // IWYU pragma: keep
#include <iostream> // IWYU pragma: keep
#define PHARE_LOG_LINE_STR(str)                                                                    \
    std::cout << __FILE__ << ":" << __LINE__ << " - " << str << std::endl;
#define PHARE_LOG_LINE_SS(s) PHARE_LOG_LINE_STR((std::stringstream{} << s).str());
#else
#define PHARE_LOG_LINE_STR(str)
#define PHARE_LOG_LINE_SS(str)
#endif
#define PHARE_LOG_LINE PHARE_LOG_LINE_STR("")


#if PHARE_WITH_CALIPER
#include "caliper/cali.h" // IWYU pragma: keep


#define PHARE_LOG_START(lvl, str) CALI_MARK_BEGIN(str)
#define PHARE_LOG_STOP(lvl, str) CALI_MARK_END(str)
#define PHARE_LOG_SCOPE(lvl, str) PHARE::scope_log PHARE_STR_CAT(__phare_scope, __LINE__)(lvl, str)

#else // !PHARE_WITH_CALIPER

#include "core/utilities/logger/logger_defaults.hpp" // IWYU pragma: keep


#endif // PHARE_WITH_CALIPER

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


#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO) || defined(PHARE_FORCE_LOG_LINE)
struct ScopeTimer
{
    std::string key;
    std::size_t start_time = now();

    std::string static inline const FORMAT = ":%Y-%m-%d-%H:%M:%S";

    ~ScopeTimer()
    {
        std::cout << formated_time() << " PHARE SCOPE TIMER: " << key
                  << " time: " << (now() - start_time) / 1e6 << " ms" << std::endl;
    }

    std::uint64_t static now()
    {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(
                   std::chrono::steady_clock::now().time_since_epoch())
            .count();
    }

    std::string static formated_time()
    {
        char date[256];
        auto now       = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        std::strftime(date, sizeof(date), FORMAT.c_str(), std::gmtime(&in_time_t));
        return date;
    }
};
#endif //


} // namespace PHARE

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO) || defined(PHARE_FORCE_LOG_LINE)
#define PHARE_FN_TIMER(key) PHARE::ScopeTimer PHARE_STR_CAT(__phare_scope_, __LINE__)(key)
#else
#define PHARE_FN_TIMER(key) // noop
#endif                      //

#endif /* PHARE_CORE_LOGGER_H */
