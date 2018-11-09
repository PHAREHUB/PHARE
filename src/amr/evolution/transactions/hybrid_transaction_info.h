
#ifndef PHARE_HYBRID_TRANSACTION_INFO_H
#define PHARE_HYBRID_TRANSACTION_INFO_H

#include "transaction_info.h"


#include <string>
#include <vector>




namespace PHARE
{
/**
 * @brief The HybridTransactionInfo class derives from ITransactionInfo. It is used to setup a
 * HybridTransaction and thus contains all the names of all variables that will need to be
 * communicated by a HybridTransaction.
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


class HybridTransactionInfo : public ITransactionInfo
{
public:
    VecFieldNames modelMagnetic;
    VecFieldNames modelElectric;

    std::vector<VecFieldNames> initMagnetic;
    std::vector<VecFieldNames> initElectric;

    std::vector<VecFieldNames> ghostMagnetic;
    std::vector<VecFieldNames> ghostElectric;

    virtual ~HybridTransactionInfo() = default;
};



} // namespace PHARE
#endif
