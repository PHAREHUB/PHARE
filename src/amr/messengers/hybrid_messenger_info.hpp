#ifndef PHARE_HYBRID_MESSENGER_INFO_HPP
#define PHARE_HYBRID_MESSENGER_INFO_HPP

#include "messenger_info.hpp"

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
    public:
        // store names of field and vector fields known to be part of the model
        // i.e. that constitute the state of the model between two time steps.
        std::string modelMagnetic;
        std::string modelElectric;
        std::string modelIonChargeFlux;
        std::string modelIonDensity;

        // fully kinetic model only: the particle current, deposited by particles and summed on
        // patch borders
        std::string modelCurrent;

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


        std::vector<std::string> maxBorderFields;
        std::vector<std::string> maxBorderVecFields;

        // below are the descriptions of the electric field that we use in the refluxing
        std::string refluxElectric;
        std::string fluxSumElectric;

        virtual ~HybridMessengerInfo() = default;
    };


    /**
     * @brief fillParticleModelMessengerInfo declares the quantities shared by models whose species
     * are particles (hybrid and fully kinetic): EM fields, particles and moments.
     */
    void fillParticleModelMessengerInfo(auto const& state, amr::HybridMessengerInfo& info)
    {
        // only the charge density and the charge flux, which give the electron velocity in hybrid,
        // are communicated. The mass density and the bulk velocity are diagnostics only
        info.modelMagnetic      = state.electromag.B.name();
        info.modelElectric      = state.electromag.E.name();
        info.modelIonDensity    = state.ions.chargeDensityName();
        info.modelIonChargeFlux = state.ions.chargeFlux().name();

        info.initElectric.emplace_back(state.electromag.E.name());
        info.initMagnetic.emplace_back(state.electromag.B.name());

        info.ghostElectric.push_back(info.modelElectric);
        info.ghostMagnetic.push_back(info.modelMagnetic);

        for (auto const& pop : state.ions)
        {
            info.interiorParticles.push_back(pop.name());
            info.levelGhostParticlesOld.push_back(pop.name());
            info.levelGhostParticlesNew.push_back(pop.name());
            info.patchGhostParticles.push_back(pop.name());
            info.ghostFlux.emplace_back(pop.flux().name());
            info.sumBorderFields.emplace_back(pop.particleDensity().name());
        }

        info.maxBorderFields.emplace_back(state.ions.chargeDensity().name());
        info.maxBorderVecFields.emplace_back(state.ions.chargeFlux().name());
    }


} // namespace amr
} // namespace PHARE
#endif
