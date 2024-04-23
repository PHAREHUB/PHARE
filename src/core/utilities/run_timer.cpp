
#include "run_timer.hpp"

namespace PHARE::core
{

namespace detail
{

    static bool _scope_timing_active = false;

} // namespace detail

StaticRunTimer& StaticRunTimer::INSTANCE()
{
    static StaticRunTimer i;
    std::string const e          = core::get_env("PHARE_SCOPE_TIMING", "false");
    detail::_scope_timing_active = e == "1" || e == "true";
    return i;
}

namespace detail
{
    static scope_timer* _current_scope_timer = nullptr;

} // namespace detail

scope_timer::scope_timer(RunTimerReport& _r)
    : r{_r}
{
    if (detail::_scope_timing_active)
    {
        this->pscope = detail::_current_scope_timer;

        if (this->pscope)
            pscope->childs.reserve(pscope->childs.size() + 1);

        detail::_current_scope_timer = this;
    }
}

scope_timer::~scope_timer()
{
    if (detail::_scope_timing_active)
    {
        detail::_current_scope_timer = this->pscope;

        auto& s = *r.snapshots.emplace_back( // allocated in construtor
            std::make_shared<RunTimerReportSnapshot>(&r, parent, now() - start));

        if (this->pscope)
            pscope->childs.emplace_back(&s);

        s.childs = std::move(childs);

        if (parent == nullptr)
            StaticRunTimer::INSTANCE().traces.emplace_back(&s);

        StaticRunTimer::INSTANCE().report_stack_ptr = parent;
    }
}


} // namespace PHARE::core
