#ifndef PHARE_CORE_ENV_HPP
#define PHARE_CORE_ENV_HPP

// Single source for handling env vars

#include <array>
#include <cstdint>
#include <optional>
#include <string_view>
#include <unordered_map>
#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"

namespace PHARE::env
{

struct Var
{
    using value_type                   = std::string;
    using results_type                 = std::unordered_map<std::string, std::string>;
    auto constexpr static noop_results = [](Var const&) { return results_type{}; };

    std::string_view const id;
    std::string_view const desc;
    std::vector<std::pair<std::string_view, std::string_view>> const options;

    std::optional<value_type> const _default               = std::nullopt;
    std::function<results_type(Var const&)> const _results = noop_results;
    std::optional<value_type> const v                      = get();
    results_type const results                             = _results(*this);

    std::optional<value_type> get() const
    {
        std::string _id{id};
        if (_default)
            return core::get_env(_id, *_default);
        return core::get_env(_id);
    }
    auto const& operator()() const { return v; }
    auto const& operator()(std::string const& s) const { return results.at(s); }
    auto const& operator()(std::string const& s, std::string const& default_) const
    {
        if (results.count(s))
            return results.at(s);
        return default_;
    }
    bool exists() const { return v != std::nullopt; }
    operator bool() const { return exists(); }
};

} // namespace PHARE::env

namespace PHARE
{
class Env
{
public:
    template<typename T>
    using map_t = std::unordered_map<std::string, T const* const>;

    static Env& INSTANCE()
    {
        if (!self)
            self = std::make_unique<Env>();
        return *self;
    }
    static auto& reinit() { return *(self = std::make_unique<Env>()); }

    env::Var const PHARE_LOG{
        "PHARE_LOG",
        "Write logs to $CWD/.log",
        {{{"RANK_FILES", "Write logs $CWD/.log, a file per rank"},
          {"DATETIME_FILES", "Write logs $CWD/.log, filename per rank and datetime"},
          {"NONE", "print normally to std::cout"}}},
        std::nullopt,
        [](auto const& self) {
            std::string static const file_key = "PHARE_LOG_FILE";
            typename env::Var::results_type map;
            if (auto const& opt = self())
            {
                auto const& val = *opt;
                if (val == "RANK_FILES")
                    map[file_key] = ".log/" + std::to_string(core::mpi::rank()) + ".out";
                else if (val == "DATETIME_FILES")
                {
                    auto date_time = core::mpi::date_time();
                    auto rank      = std::to_string(core::mpi::rank());
                    auto size      = std::to_string(core::mpi::size());
                    map[file_key]  = ".log/" + date_time + "_" + rank + "_of_" + size + ".out";
                }
                else if (val != "NONE")
                    throw std::runtime_error("PHARE_LOG invalid type, valid keys are "
                                             "RANK_FILES/DATETIME_FILES/NONE");
            }
            return map;
        } //
    };
    env::Var const PHARE_SCOPE_TIMING{
        "PHARE_SCOPE_TIMING", "Enable function scope timing", {{{"1", "ON"}, {"0", "OFF"}}}, "0"};

    map_t<env::Var> const vars = {{
        {"PHARE_LOG", &PHARE_LOG},
        {"PHARE_SCOPE_TIMING", &PHARE_SCOPE_TIMING},
    }};

    auto& operator()(std::string const& s) const
    {
        assert(vars.count(s));
        return *vars.at(s);
    }

private:
    static inline std::unique_ptr<Env> self = nullptr;
};

} // namespace PHARE

#endif /* PHARE_CORE_ERRORS_H */
