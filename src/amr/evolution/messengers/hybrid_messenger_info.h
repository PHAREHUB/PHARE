
#ifndef PHARE_HYBRID_MESSENGER_INFO_H
#define PHARE_HYBRID_MESSENGER_INFO_H

#include "messenger_info.h"


#include <string>
#include <vector>




namespace PHARE
{
/**
 * @brief The HybridMessengerInfo class derives from IMessengerInfo. It is used to setup a
 * HybridMessenger and thus contains all the names of all variables that will need to be
 * communicated by a HybridMessenger.
 *
 * One can distinguish two groups of variables:
 *
 * - init variables, are the names of the variables that need to be initialized
 * - ghost variables are the names of the variables that need to be filled in ghost zones
 *
 * these variables can belong to the model, the solver.
 *
 */

struct VecFieldNames
{
    std::string vecName;
    std::string xName;
    std::string yName;
    std::string zName;
};


class HybridMessengerInfo : public IMessengerInfo
{
public:
    VecFieldNames modelMagnetic;
    VecFieldNames modelElectric;

    std::vector<VecFieldNames> initMagnetic;
    std::vector<VecFieldNames> initElectric;

    std::vector<VecFieldNames> ghostMagnetic;
    std::vector<VecFieldNames> ghostElectric;

    virtual ~HybridMessengerInfo() = default;
};



} // namespace PHARE
#endif
