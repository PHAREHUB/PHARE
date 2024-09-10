
#ifndef PHARE_MHD_MESSENGER_INFO_HPP
#define PHARE_MHD_MESSENGER_INFO_HPP

#include "core/data/vecfield/vecfield.hpp"
#include "messenger_info.hpp"



namespace PHARE
{
namespace amr
{
    class MHDMessengerInfo : public IMessengerInfo
    {
        using VecFieldNames = core::VecFieldNames;
    public:
        // What is needed here ?
        virtual ~MHDMessengerInfo() = default;
    };

} // namespace amr


} // namespace PHARE
#endif
