#ifndef PHARE_CORE_UTILITIES_META_META_UTILITIES_H
#define PHARE_CORE_UTILITIES_META_META_UTILITIES_H

#include <iterator>
#include <type_traits>

namespace PHARE
{
template<typename...>
using tryToInstanciate = void;


struct dummy
{
    using type              = int;
    static const type value = 0;
};




template<typename IterableCandidate, typename AttemptBegin = void, typename AttemptEnd = void>
struct has_beginend : std::false_type
{
};



/** \brief has_beginend is a traits that permit to check if a Box or a BoxContainer
 * is passed as template argument
 */
template<typename IterableCandidate>
struct has_beginend<IterableCandidate,
                    tryToInstanciate<decltype(std::begin(std::declval<IterableCandidate>()))>,
                    tryToInstanciate<decltype(std::end(std::declval<IterableCandidate>()))>>
    : std::true_type
{
};


template<typename IterableCandidate>
using is_iterable = std::enable_if_t<has_beginend<IterableCandidate>::value, dummy::type>;

} // namespace PHARE




#endif
