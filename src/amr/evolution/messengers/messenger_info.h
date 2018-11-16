
#ifndef PHARE_MESSENGER_INFO_H
#define PHARE_MESSENGER_INFO_H

#include <string>



namespace PHARE
{
/**
 * @brief The IMessengerInfo class is an abstract class used to setup an abstract IMessenger
 * object. This class needs to be abstract since an IMessenger does not know which concrete kind
 * of messenger information is needed. Derived implemetation should provide all informations
 * required by a concrete implementation of IMessenger::setup()
 */
class IMessengerInfo
{
public:
    virtual ~IMessengerInfo() = default;
};




} // namespace PHARE
#endif
