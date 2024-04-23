#ifndef PHARE_CORE_UTILITIES_RUN_TIMER_HPP
#define PHARE_CORE_UTILITIES_RUN_TIMER_HPP

#include <array>
#include <memory>
#include <vector>
#include <thread>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string_view>
#include <unordered_map>

#include "core/logger.hpp"
#include "core/utilities/mpi_utils.hpp"


namespace PHARE::core
{

struct RunTimerReport;
struct RunTimerReportSnapshot;

namespace detail
{

    void inline write_timer_file();

} // namespace detail

struct StaticRunTimer
{
    static StaticRunTimer& INSTANCE();

    StaticRunTimer()
        : _rank{mpi::is_init() ? mpi::rank() : 0}
        , timer_file{".phare_times." + std::to_string(_rank) + ".bin"}
    {
        reports.reserve(5);
        traces.reserve(5);
    }

    void init() {}

    void shutdown()
    {
        if (traces.size())
        {
            detail::write_timer_file();
            traces.clear();
            reports.clear();
        }
    }

    static void reset() { INSTANCE().shutdown(); }

    std::vector<RunTimerReport*> reports;
    std::vector<RunTimerReportSnapshot*> traces;
    RunTimerReport* report_stack_ptr = nullptr;
    std::uint32_t n_reports          = 0;

    int const _rank = 0;
    std::string const timer_file;
};


struct RunTimerReportSnapshot
{
    RunTimerReportSnapshot(RunTimerReport* s, RunTimerReport* p, std::uint32_t&& t)
        : self{s}
        , parent{p}
        , time{t}
    {
        childs.reserve(2);
    }

    RunTimerReport const* const self;
    RunTimerReport const* const parent;

    std::uint64_t const time;
    std::vector<RunTimerReportSnapshot*> childs;
};

struct RunTimerReport
{
    std::string_view k, f;
    std::uint32_t l = 0;

    RunTimerReport(std::string_view const& _k, std::string_view const& _f, std::uint32_t const& _l)
        : k{_k}
        , f{_f}
        , l{_l}
    {
        StaticRunTimer::INSTANCE().reports.emplace_back(this);
        ++StaticRunTimer::INSTANCE().n_reports;
        snapshots.reserve(5);
    }

    ~RunTimerReport() {}

    auto operator()(std::size_t i) { return snapshots[i].get(); }
    auto size() { return snapshots.size(); }

    std::vector<std::shared_ptr<RunTimerReportSnapshot>> snapshots; // emplace back breaks pointers!
};



static void inline run_timer_shutdown()
{
    StaticRunTimer::INSTANCE().shutdown();
}



struct scope_timer
{
    scope_timer(RunTimerReport& _r);

    ~scope_timer();

    std::uint64_t static now()
    {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(
                   std::chrono::system_clock::now().time_since_epoch())
            .count();
    }

    static scope_timer& root_parent_from(scope_timer& self)
    {
        if (self.pscope)
            return root_parent_from(*self.pscope);
        return self;
    }


    RunTimerReport& r;
    RunTimerReport* parent = StaticRunTimer::INSTANCE().report_stack_ptr;
    RunTimerReport* child  = nullptr;

    scope_timer* pscope       = nullptr;
    std::uint64_t const start = now();

    std::vector<RunTimerReportSnapshot*> childs;
};


struct BinaryTimerFileNode
{
    BinaryTimerFileNode(std::uint16_t _fn_id, std::uint64_t _time)
        : fn_id{_fn_id}
        , time{_time}
    {
    }

    std::uint16_t fn_id;
    std::uint64_t time;

    std::vector<BinaryTimerFileNode> kinder{};
};

struct BinaryTimerFile
{
    BinaryTimerFile(bool fromStaticTraces = true)
    {
        if (fromStaticTraces)
        {
            for (auto const& trace : StaticRunTimer::INSTANCE().traces)
                recurse_traces_for_keys(trace);
            for (auto const& trace : StaticRunTimer::INSTANCE().traces)
                recurse_traces_for_nodes(
                    trace, roots.emplace_back(key_ids[std::string{trace->self->k}], trace->time));
        }
    }


    template<typename Trace>
    void recurse_traces_for_nodes(Trace const& c, BinaryTimerFileNode& node)
    {
        for (std::size_t i = 0; i < c->childs.size(); ++i)
            recurse_traces_for_nodes(
                c->childs[i], node.kinder.emplace_back(key_ids[std::string{c->childs[i]->self->k}],
                                                       c->childs[i]->time));
    }

    template<typename Trace>
    void recurse_traces_for_keys(Trace const& c)
    {
        std::string s{c->self->k};
        if (!key_ids.count(s))
        {
            auto [it, b] = key_ids.emplace(s, key_ids.size());
            assert(b);
            auto const& [k, i] = *it;
            assert(!id_to_key.count(i));
            id_to_key.emplace(i, k);
        }
        for (std::size_t i = 0; i < c->childs.size(); ++i)
            recurse_traces_for_keys(c->childs[i]);
    }

    template<typename... Args>
    void _byte_write(std::ofstream& file, Args const... args) const
    {
        std::stringstream ss;
        (ss << ... << args);
        auto s = ss.str();

        file.write(s.c_str(), s.size());
        file << std::endl;
    }

    void write(std::string const& filename) const
    {
        std::ofstream f{filename, std::ios::binary};
        for (auto const& [i, k] : id_to_key)
            _byte_write(f, i, " ", k);
        f << std::endl; // break between function id map and function times
        for (auto const& root : roots)
            _write(f, root);
    }

    void _write(std::ofstream& file, BinaryTimerFileNode const& node, std::uint16_t tabs = 0) const
    {
        for (std::size_t ti = 0; ti < tabs; ++ti)
            file << " ";
        file << node.fn_id << " " << node.time << std::endl;
        for (auto const& n : node.kinder)
            _write(file, n, tabs + 1);
    }


    std::unordered_map<std::string, std::size_t> key_ids; // only used on write
    std::unordered_map<std::size_t, std::string> id_to_key;
    std::vector<BinaryTimerFileNode> roots;
};


namespace detail
{

    void inline write_timer_file()
    {
        BinaryTimerFile{}.write(StaticRunTimer::INSTANCE().timer_file);
    }


} // namespace detail


} // namespace PHARE::core


#define RUN_TIMER_SCOPE(key)                                                                       \
    static PHARE::core::RunTimerReport STR_CAT(ridx_, __LINE__){key, __FILE__, __LINE__};          \
    PHARE::core::scope_timer STR_CAT(_scope_timer_, __LINE__){STR_CAT(ridx_, __LINE__)};           \
    PHARE::core::StaticRunTimer::INSTANCE().report_stack_ptr = &STR_CAT(ridx_, __LINE__);


#endif /*PHARE_CORE_UTILITIES_RUN_TIMER_HPP*/
