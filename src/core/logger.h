

#include <tuple>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>


#include "thread_pool.h"

auto get_thread_id()
{
    std::ostringstream os;
    os << std::hex << pthread_self();
    return os.str();
}

auto now_in_nanos()
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
        pool.enqueue([this](std::string s1) { this->writer << s1 << std::endl; }, s0);
    }
    void log(std::string&& s) { log(s); }

    ThreadPool pool{1}; // overkill but useful for the moment
    std::ofstream writer{std::string{"Logger."} + get_thread_id() + ".txt"};
    std::vector<StartStopLogger> nestings;
};

class LogMan
{
public:
    LogMan(std::string base_path)
        : base_path_{base_path}
    {
    }

    auto _fn() { return std::forward_as_tuple(get_thread_id()); }


    void start(std::string&& key) { loggers_.at(get_thread_id())->start(key); }
    void stop() { loggers_.at(get_thread_id())->stop(); }
    auto scope(std::string&& key) { return loggers_.at(get_thread_id())->scope(key); }


    LogMan& registerLogger()
    {
        auto const thread_id = get_thread_id();
        if (!loggers_.count(thread_id))
            loggers_.emplace(thread_id, std::make_shared<Logger>());
        return *this;
    }

    // void log(std::string const& s) { loggers_.at(get_thread_id())->log(s); }

private:
    std::string base_path_;
    std::unordered_map<std::string /*thread_id*/, std::shared_ptr<Logger>> loggers_;
};
