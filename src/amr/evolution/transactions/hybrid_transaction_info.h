
#ifndef PHARE_HYBRID_TRANSACTION_INFO_H
#define PHARE_HYBRID_TRANSACTION_INFO_H

#include "transaction_info.h"


#include <string>
#include <vector>




namespace PHARE
{
// une hybride transaction info est donnée par un modele hybride
// ou un solver hybride pour que la transaction créé les algos dont elle a besoin
//
// ce qui en fait une Hybridxxxx c'est qu'elle sait ce que sont que
// des champs et des particules etc.

// elle va créer les refineAlgo donc elle a besoin de connaitre
// si temporel ou non temporel
//      si non temporel:
//          - variable source et de destination
//        - variables old et new de depart et variable de destination

class HybridTransactionInfo : public ITransactionInfo
{
public:
    std::string modelMagneticX;
    std::string modelMagneticY;
    std::string modelMagneticZ;
    std::string modelMagneticName;

    std::string modelElectricX;
    std::string modelElectricY;
    std::string modelElectricZ;
    std::string modelElectricName;

    std::vector<std::string> solverMagneticX;
    std::vector<std::string> solverMagneticY;
    std::vector<std::string> solverMagneticZ;
    std::vector<std::string> solverMagneticName;

    std::vector<std::string> solverElectricX;
    std::vector<std::string> solverElectricY;
    std::vector<std::string> solverElectricZ;
    std::vector<std::string> solverElectricName;




    // we don't need the model electric, since the electric field is assign
    // by the solver anyway

    virtual ~HybridTransactionInfo() = default;
};



} // namespace PHARE
#endif
