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

    void static LOG(std::string s)
    {
        for (std::string const needle : {"PHARE::core::", "PHARE::"})
            for (auto pos = s.find(needle); pos != std::string::npos; pos = s.find(needle))
                s.erase(pos, needle.size());
        INSTANCE().ops[s] += 1;
    }

    std::size_t static COUNT(std::string const& s)
    {
        auto const& ops = INSTANCE().ops;
        auto const it   = ops.find(s);
        return it == ops.end() ? 0 : it->second;
    }

    void static RESET() { INSTANCE().ops.clear(); }

    auto static MOVE()
    {
        MemoryMonitoring mm{std::move(INSTANCE().ops)};
        INSTANCE().ops = {};
        return mm;
    }

    void print() const
    {
        for ([[maybe_unused]] auto const& [k, v] : ops)
        {
            PHARE_LOG_LINE_SS(k << " " << v);
        }
    }

    void static PRINT() { INSTANCE().print(); }

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
