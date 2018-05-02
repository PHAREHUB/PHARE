#ifndef PHARE_CORE_UTILITIES_META_META_UTILITIES_H
#define PHARE_CORE_UTILITIES_META_META_UTILITIES_H

namespace PHARE
{
template<typename...>
using tryToInstanciate = void;
}

struct dummy
{
    using type              = int;
    static const type value = 0;
};

#endif
