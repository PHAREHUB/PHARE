#include "messenger_factory.hpp"

namespace PHARE::amr
{
std::vector<MessengerDescriptor> makeDescriptors(std::vector<std::string> modelNames)
{
    if (modelNames.size() == 1)
        return {{modelNames[0], modelNames[0]}};

    else if (modelNames.size() == 2)
    {
        return {{modelNames[0], modelNames[0]},
                {modelNames[0], modelNames[1]},
                {modelNames[1], modelNames[1]}};
    }
    else
        throw std::runtime_error("Error max number of models is 2");
}
} // namespace PHARE::amr
