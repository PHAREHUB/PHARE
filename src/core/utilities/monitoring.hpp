#ifndef PHARE_CORE_UTILITIES_MONITORING_HPP
#define PHARE_CORE_UTILITIES_MONITORING_HPP

#include "core/logger.hpp"
#include "core/def/phare_config.hpp"

#include <string>
#include <cstddef>
#include <unordered_map>


namespace PHARE::core
{

// not thread safe!
struct MemoryMonitoring
{
    static auto& INSTANCE()
    {
        static MemoryMonitoring i;
        return i;
    }

    void static LOG(std::string const& s) { INSTANCE().ops[s] += 1; }

    void static PRINT()
    {
        for ([[maybe_unused]] auto const& [k, v] : INSTANCE().ops)
        {
            PHARE_LOG_LINE_SS(k << " " << v);
        }
    }

    std::unordered_map<std::string, std::size_t> ops;
};



struct MemoryMonitor
{
    void create() { MemoryMonitoring::LOG(s + "::construct"); }
    void copy() { MemoryMonitoring::LOG(s + "::copy_construct"); }
    void move() { MemoryMonitoring::LOG(s + "::move_construct"); }
    void move_assign() { MemoryMonitoring::LOG(s + "::move_assign"); }
    void copy_assign() { MemoryMonitoring::LOG(s + "::copy_assign"); }

    std::string s;
};



} // namespace PHARE::core




#endif /* PHARE_CORE_UTILITIES_MONITORING_HPP */
