#include "debugod.hpp"


namespace PHARE::amr
{
template<typename PHARE_TYPES>
DEBUGOD<PHARE_TYPES>& DEBUGOD<PHARE_TYPES>::INSTANCE()
{
    static DEBUGOD instance;
    return instance;
}
}; // namespace PHARE::amr
