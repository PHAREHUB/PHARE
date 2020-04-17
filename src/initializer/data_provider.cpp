
#include "initializer/data_provider.h"


namespace PHARE::initializer
{
PHAREDictHandler& PHAREDictHandler::INSTANCE()
{
    static PHAREDictHandler handler;
    return handler;
}
} // namespace PHARE::initializer
