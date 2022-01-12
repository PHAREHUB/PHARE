#ifndef PHARE_RESTART_DATA_PROVIDER_HPP
#define PHARE_RESTART_DATA_PROVIDER_HPP

#include "initializer/data_provider.hpp"



namespace PHARE
{
namespace initializer
{
    class RestartDataProvider : public DataProvider
    {
    public:
        virtual void read() override
        {
            // read in python script
        }
    };

} // namespace initializer

} // namespace PHARE

#endif // DATA_PROVIDER_HPP
