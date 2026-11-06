#ifndef PHARE_CORE_UTILITIES_TUPLE_HPP
#define PHARE_CORE_UTILITIES_TUPLE_HPP

#include <tuple>

namespace PHARE::core
{
template<typename Tuple, typename Func>
void apply(Tuple tuple, Func func)
{
    std::apply([&](auto&... args) { (func(args), ...); }, tuple);
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_TUPLE_HPP
