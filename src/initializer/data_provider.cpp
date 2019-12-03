
#include "initializer/data_provider.h"


namespace PHARE
{
namespace initializer
{
    PHAREDictHandler& PHAREDictHandler::INSTANCE()
    {
        static PHAREDictHandler handler;
        return handler;
    }
} // namespace initializer

} // namespace PHARE
