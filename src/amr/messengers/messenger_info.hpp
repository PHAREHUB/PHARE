
#ifndef PHARE_MESSENGER_INFO_HPP
#define PHARE_MESSENGER_INFO_HPP

#include <string>



namespace PHARE
{
namespace amr
{
    /**
     * @brief The IMessengerInfo class is an abstract class used to register quantities to be used
     * by a Messenger in communications. This class is abstract and each concrete IMessenger will
     * need to define its concrete specialization of the IMessengerInfo.
     *
     * * This class and its subclasses are intermediate objects between IPhysicalModel and ISolver
     * subclasses and IMessenger subclasses.
     *
     * Concrete implementations of this class are used by concrete implementations of
     * IPhysicalModel, ISolver to declare the quantities they need concrete IMessenger to deal with.
     *
     */
    class IMessengerInfo
    {
    public:
        virtual ~IMessengerInfo() = default;
    };


} // namespace amr


} // namespace PHARE
#endif
