
#include "python_data_provider.h"


namespace PHARE
{
namespace initializer
{
    PythonDictHandler& PythonDictHandler::INSTANCE()
    {
        static PythonDictHandler handler;
        return handler;
    }

} /*namespace initializer*/
} /*namespace PHARE*/
