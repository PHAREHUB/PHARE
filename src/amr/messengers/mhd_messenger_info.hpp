#ifndef PHARE_MHD_MESSENGER_INFO_HPP
#define PHARE_MHD_MESSENGER_INFO_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "messenger_info.hpp"



namespace PHARE
{
namespace amr
{
    class MHDMessengerInfo : public IMessengerInfo
    {
    public:
        std::string modelDensity;
        std::string modelVelocity;
        std::string modelMagnetic;
        std::string modelPressure;

        std::string modelMomentum;
        std::string modelTotalEnergy;

        std::string modelElectric;
        std::string modelCurrent;

        std::vector<std::string> initDensity;
        std::vector<std::string> initMomentum;
        std::vector<std::string> initMagnetic;
        std::vector<std::string> initTotalEnergy;

        std::vector<std::string> ghostDensity;
        std::vector<std::string> ghostVelocity;
        std::vector<std::string> ghostMagnetic; // not actually to fill ghost cells but rather for
                                                // amr operations, see hybrid
        std::vector<std::string> ghostPressure;
        std::vector<std::string> ghostMomentum;
        std::vector<std::string> ghostTotalEnergy;
        std::vector<std::string> ghostMagneticFluxesX;
        std::vector<std::string> ghostMagneticFluxesY;
        std::vector<std::string> ghostMagneticFluxesZ;
        std::vector<std::string> ghostElectric;
        std::vector<std::string> ghostCurrent;

        core::AllFluxesNames reflux;
        core::AllFluxesNames fluxSum;
        std::string refluxElectric;
        std::string fluxSumElectric;

        virtual ~MHDMessengerInfo() = default;
    };

} // namespace amr


} // namespace PHARE
#endif
