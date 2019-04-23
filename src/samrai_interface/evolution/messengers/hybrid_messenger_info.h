
#ifndef PHARE_HYBRID_MESSENGER_INFO_H
#define PHARE_HYBRID_MESSENGER_INFO_H

#include "data/vecfield/vecfield_component.h"
#include "messenger_info.h"

#include <string>
#include <vector>




namespace PHARE
{
namespace amr_interface
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

    struct VecFieldDescriptor
    {
        std::string vecName;
        std::string xName;
        std::string yName;
        std::string zName;

        VecFieldDescriptor() = default;

        template<typename VecFieldT>
        VecFieldDescriptor(VecFieldT const& v)
            : vecName{v.name()}
            , xName{v.getComponentName(core::Component::X)}
            , yName{v.getComponentName(core::Component::Y)}
            , zName{v.getComponentName(core::Component::Z)}

        {
        }
    };



    using FieldDescriptor      = std::string;
    using PopulationDescriptor = std::string;


    /// template<typename VecFieldT>
    class HybridMessengerInfo : public IMessengerInfo
    {
    public:
        VecFieldDescriptor modelMagnetic;
        VecFieldDescriptor modelElectric;


        //! names of the magnetic quantities that will be communicated by
        //! HybridMessenger::initLevel() and HybridMessenger::regrid()
        std::vector<VecFieldDescriptor> initMagnetic;


        //! names of the electric quantities that will be communicated by
        //! HybridMessenger::initLevel() and HybridMessenger::regrid()
        std::vector<VecFieldDescriptor> initElectric;

        std::vector<PopulationDescriptor> interiorParticles;
        std::vector<PopulationDescriptor> levelGhostParticlesNew;
        std::vector<PopulationDescriptor> levelGhostParticlesOld;
        std::vector<PopulationDescriptor> patchGhostParticles;


        //! name of the magnetic quantities that will be communicated by
        //! HybridMessenger::fillMagneticGhosts()
        std::vector<VecFieldDescriptor> ghostMagnetic;
        // std::vector<VecFieldT*> ghostsMagnetic;


        //! name of the electric quantities that will be communicated by the
        //! HybridMessenger::fillGhostElectric
        std::vector<VecFieldDescriptor> ghostElectric;



        virtual ~HybridMessengerInfo() = default;
    };


} // namespace amr_interface
} // namespace PHARE
#endif
