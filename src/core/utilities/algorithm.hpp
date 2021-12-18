#ifndef PHARE_ALGORITHM_H
#define PHARE_ALGORITHM_H

#include "core/utilities/types.hpp"



#include <algorithm>
#include <string>

namespace PHARE
{
namespace core
{
    template<std::uint32_t lhs, std::uint32_t rhs>
    constexpr std::uint32_t max()
    {
        if constexpr (lhs < rhs)
        {
            return rhs;
        }
        else if constexpr (lhs >= rhs)
        {
            return lhs;
        }
    }


    template<typename T>
    std::string to_str(T&& t)
    {
        return t.to_str();
    }



    template<typename Container, typename ContainedT = typename Container::value_type>
    bool notIn(ContainedT& obj, Container& list)
    {
        auto sameItem = std::find_if(std::begin(list), std::end(list), [&obj](auto& currentItem) {
            return obj->name() == currentItem->name();
        });

        return sameItem == std::end(list);
    }



} // namespace core
} // namespace PHARE

#endif
