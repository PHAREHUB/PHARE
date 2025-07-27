#ifndef PHARE_HYBRID_MESSENGER_INFO_HPP
#define PHARE_HYBRID_MESSENGER_INFO_HPP

#include "messenger_info.hpp"
#include "core/data/vecfield/vecfield.hpp"

#include <string>
#include <vector>




namespace PHARE
{
namespace amr
{
    /**
     * @brief The HybridMessengerInfo class derives from IMessengerInfo.
     * It is filled by HybridModel and a Hybrid ISolver to register the quantities they need to be
     * communicated by the HybridMessenger.
     *
     * One can distinguish two groups of variables:
     *
     * - init variables, are the names of the variables that need to be initialized. Those are the
     * quantities that will be communicated during IMessenger::regrid() and IMessenger::initLevel()
     *
     * - ghost variables are the names of the variables that need to be filled in ghost zones.
     *
     * these variables can belong to the model, the solver.
     *
     */




    class HybridMessengerInfo : public IMessengerInfo
    {
        // using std::string = core::std::string;

    public:
        // store names of field and vector fields known to be part of the model
        // i.e. that constitute the state of the model between two time steps.
        std::string modelMagnetic;
        std::string modelElectric;
        std::string modelCurrent;
        std::string modelIonBulkVelocity;
        std::string modelIonDensity;

        // store names of vector fields that need to be initialized by refinement
        // moments are initialized by particles so only EM fields need to be init.
        std::vector<std::string> initMagnetic;
        std::vector<std::string> initElectric;

        // below are the names of the populations that need to be communicated
        // this is for initialization
        std::vector<std::string> interiorParticles;

        // and the following are for ghost cells
        std::vector<std::string> levelGhostParticlesNew;
        std::vector<std::string> levelGhostParticlesOld;
        std::vector<std::string> patchGhostParticles;

        // below are the descriptions of the vector fields that for which
        // ghosts need to be filled at some point.
        std::vector<std::string> ghostFlux;
        std::vector<std::string> sumBorderFields;
        std::vector<std::string> ghostMagnetic;
        std::vector<std::string> ghostElectric;
        std::vector<std::string> ghostCurrent;
        std::vector<std::string> ghostBulkVelocity;

        virtual ~HybridMessengerInfo() = default;
    };


} // namespace amr
} // namespace PHARE
#endif
