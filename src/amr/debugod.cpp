#include "debugod.hpp"


namespace PHARE::amr
{
DEBUGOD& DEBUGOD::INSTANCE()
{
    static DEBUGOD instance;
    return instance;
}
}; // namespace PHARE::amr
