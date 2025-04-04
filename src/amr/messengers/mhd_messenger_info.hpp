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
        std::vector<VecFieldNames> initMomentum;
        std::vector<VecFieldNames> initMagnetic;
        std::vector<std::string> initTotalEnergy;

        std::vector<std::string> ghostDensity;
        std::vector<VecFieldNames> ghostVelocity;
        std::vector<VecFieldNames> ghostMagnetic; // not actually to fill ghost cells but rather for
                                                  // amr operations, see hybrid
        std::vector<std::string> ghostPressure;
        std::vector<VecFieldNames> ghostMomentum;
        std::vector<std::string> ghostTotalEnergy;
        std::vector<VecFieldNames> ghostMagneticFluxesX;
        std::vector<VecFieldNames> ghostMagneticFluxesY;
        std::vector<VecFieldNames> ghostMagneticFluxesZ;
        std::vector<VecFieldNames> ghostElectric;
        std::vector<VecFieldNames> ghostCurrent;

        virtual ~MHDMessengerInfo() = default;
    };

} // namespace amr


} // namespace PHARE
#endif
