
#ifndef PHARE_CORE_LOGGER_H
#define PHARE_CORE_LOGGER_H


#include <tuple>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>


#include "thread_queue.h"
#include "initializer/data_provider.h"

namespace PHARE::core
{
inline auto get_thread_id()
{
    std::ostringstream os;
    os << std::hex << pthread_self();
    return os.str();
}

inline auto now_in_nanos()
{
    return static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        std::chrono::system_clock::now().time_since_epoch())
                                        .count());
}

class Logger
{
public:
    Logger() = default;

    struct StartStopLogger
    {
        StartStopLogger(std::string const& key)
            : key_{key}
        {
        }
        std::string key_;
        std::size_t start_{now_in_nanos()};
    };

    void start(std::string const& key) { nestings.emplace_back(key); }

    void stop()
    {
        auto& startStop = nestings.back();
        log(indent(nestings.size()) + startStop.key_ + " " + Logger::mtimer(startStop.start_));
        nestings.pop_back();
    }

    struct ScopeLogger
    {
        ScopeLogger(Logger& outer, std::string const& key)
            : outer_{outer}
        {
            outer_.start(key);
        }
        ~ScopeLogger() { outer_.stop(); }

        Logger& outer_;
    };

    auto scope(std::string const& key) { return ScopeLogger{*this, key}; }

private:
    Logger(Logger const&) = delete;
    Logger(Logger&&)      = delete;

    static std::string mtimer(std::size_t nanos) { return std::to_string(now_in_nanos() - nanos); }
    static std::string indent(std::size_t size) { return std::string(size, ' '); }

    void log(std::string const& s0)
    {
        queue.enqueue([this](std::string s1) { this->writer << s1 << std::endl; }, s0);
    }
    void log(std::string&& s) { log(s); }

    ThreadQueue queue; // overkill but useful for the moment
    std::ofstream writer{std::string{"Logger."} + get_thread_id() + ".txt"};
    std::vector<StartStopLogger> nestings;
};

class LogMan
{
public: // static
    static LogMan& start()
    {
        if (!self)
            self = std::make_unique<LogMan>("data");
        return *self;
    }
    static LogMan& start(PHARE::initializer::PHAREDict const& dict)
    {
        if (!self)
            self = std::make_unique<LogMan>("data");
        return *self;
    }
    static void kill() { self.release(); }

    static LogMan& get() { return *self; }

public:
    LogMan(std::string base_path)
        : base_path_{base_path}
    {
    }

    void start(std::string&& key) { logger.start(key); }
    void stop() { logger.stop(); }
    auto scope(std::string&& key) { return logger.scope(key); }

private:
    Logger logger;
    std::string base_path_;

    static std::unique_ptr<LogMan> self;
};

} // namespace PHARE::core

#define PHARE_LOG_STR() std::string{__FILE__} + " " + std::to_string(__LINE__) + " "
#define PHARE_LOG_START(str) PHARE::core::LogMan::get().start(PHARE_LOG_STR() + str)
#define PHARE_LOG_STOP() PHARE::core::LogMan::get().stop()
#define PHARE_LOG_SCOPE(str)                                                                       \
    auto _scopeLog = PHARE::core::LogMan::get().scope(PHARE_LOG_STR() + str)
#define PHARE_LOG_NAMED_SCOPE(name, str)                                                           \
    auto name = PHARE::core::LogMan::get().scope(PHARE_LOG_STR() + str)


#endif /* PHARE_CORE_LOGGER_H */
