#ifndef PHARE_RESTART_DATA_PROVIDER_H
#define PHARE_RESTART_DATA_PROVIDER_H

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

#endif // DATA_PROVIDER_H
