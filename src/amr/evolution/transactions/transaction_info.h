
#ifndef PHARE_TRANSACTION_INFO_H
#define PHARE_TRANSACTION_INFO_H

#include <string>



namespace PHARE
{
/**
 * @brief The ITransactionInfo class is an abstract class used to setup an abstract ITransaction
 * object. This class needs to be abstract since an ITransaction does not know which concrete kind
 * of transaction information is needed. Derived implemetation should provide all informations
 * required by a concrete implementation of ITransaction::setup()
 */
class ITransactionInfo
{
public:
    virtual ~ITransactionInfo() = default;
};




} // namespace PHARE
#endif
