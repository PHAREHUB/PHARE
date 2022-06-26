
#include "initializer/data_provider.hpp"


namespace PHARE::initializer
{
PHAREDictHandler& PHAREDictHandler::INSTANCE()
{
    static PHAREDictHandler handler;
    return handler;
}
} // namespace PHARE::initializer
