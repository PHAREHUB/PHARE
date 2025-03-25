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
        std::string modelDensity;
        VecFieldNames modelVelocity;
        VecFieldNames modelMagnetic;
        std::string modelPressure;

        VecFieldNames modelMomentum;
        std::string modelTotalEnergy;

        VecFieldNames modelElectric;
        VecFieldNames modelCurrent;

        std::vector<std::string> initDensity;
        std::vector<VecFieldNames> initVelocity;
        std::vector<VecFieldNames> initMagnetic;
        std::vector<std::string> initPressure;

        std::vector<std::string> ghostDensity;
        std::vector<VecFieldNames> ghostVelocity;
        std::vector<std::string> ghostPressure;
        std::vector<VecFieldNames> ghostMagneticFluxesX;
        std::vector<VecFieldNames> ghostMagneticFluxesY;
        std::vector<VecFieldNames> ghostMagneticFluxesZ;
        std::vector<VecFieldNames> ghostElectric;

        virtual ~MHDMessengerInfo() = default;
    };

} // namespace amr


} // namespace PHARE
#endif
