
#ifndef PHARE_MHD_MESSENGER_INFO_H
#define PHARE_MHD_MESSENGER_INFO_H

#include "messenger_info.h"



namespace PHARE
{
namespace amr_interface
{
    class MHDMessengerInfo : public IMessengerInfo
    {
    public:
        virtual ~MHDMessengerInfo() = default;
    };

} // namespace amr_interface


} // namespace PHARE
#endif
