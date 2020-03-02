
#ifndef PHARE_SIMULATOR_SIMULATOR_H
#define PHARE_SIMULATOR_SIMULATOR_H

#include "include.h"
#include "phare_types.h"
#include "amr/wrappers/hierarchy.h"
#include "amr/wrappers/integrator.h"



namespace PHARE
{
class ISimulator
{
public:
    virtual void initialize()       = 0;
    virtual double startTime()      = 0;
    virtual double endTime()        = 0;
    virtual double timeStep()       = 0;
    virtual void advance(double dt) = 0;
    virtual std::string to_str()    = 0;
    virtual ~ISimulator() {}
};




template<std::size_t _dimension, std::size_t _interp_order>
class Simulator : public ISimulator
{
public:
    static constexpr size_t dimension    = _dimension;
    static constexpr size_t interp_order = _interp_order;

    using SAMRAITypes = PHARE::amr::SAMRAI_Types;
    using PHARETypes  = PHARE_Types<dimension, interp_order>;

    using IPhysicalModel = PHARE::solver::IPhysicalModel<SAMRAITypes>;
    using HybridModel    = typename PHARETypes::HybridModel_t;
    using MHDModel       = typename PHARETypes::MHDModel_t;

    using SolverMHD = typename PHARETypes::SolverMHD_t;
    using SolverPPC = typename PHARETypes::SolverPPC_t;

    using MessengerFactory       = typename PHARETypes::MessengerFactory;
    using MultiPhysicsIntegrator = typename PHARETypes::MultiPhysicsIntegrator;


    Simulator(PHARE::initializer::PHAREDict dict,
              std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy)
        : hierarchy_{hierarchy}
        , modelNames_{"HybridModel"}
        , descriptors_{PHARE::amr::makeDescriptors(modelNames_)}
        , messengerFactory_{descriptors_}
        , maxLevelNumber_{dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>()}
        , dt_{dict["simulation"]["time_step"].template to<double>()}
        , timeStepNbr_{dict["simulation"]["time_step_nbr"].template to<int>()}
        , finalTime_{dt_ * timeStepNbr_}
        , multiphysInteg_{std::make_shared<MultiPhysicsIntegrator>(
              dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>())}
    {
        if (find_model("HybridModel"))
        {
            hybridModel_ = std::make_shared<HybridModel>(
                dict["simulation"],
                std::make_shared<typename HybridModel::resources_manager_type>());


            hybridModel_->resourcesManager->registerResources(hybridModel_->state);

            // we register the hybrid model for all possible levels in the hierarchy
            // since for now it is the only model available
            // same for the solver
            multiphysInteg_->registerModel(0, maxLevelNumber_ - 1, hybridModel_);
            multiphysInteg_->registerAndInitSolver(
                0, maxLevelNumber_ - 1, std::make_unique<SolverPPC>(dict["simulation"]["algo"]));
            multiphysInteg_->registerAndSetupMessengers(messengerFactory_);


            auto startTime = 0.; // TODO make it runtime
            auto endTime   = 0.; // TODO make it runtime


            integrator_ = std::make_unique<PHARE::amr::Integrator>(
                dict, hierarchy, multiphysInteg_, multiphysInteg_, startTime, endTime);
        }
        else
            throw std::runtime_error("unsupported model");
    }




    std::string to_str() override
    {
        std::stringstream ss;
        ss << "PHARE SIMULATOR\n";
        ss << "------------------------------------\n";
        ss << "interpolation order  : " << interp_order << "\n";
        ss << "dimension            : " << dimension << "\n";
        ss << "time step            : " << dt_ << "\n";
        ss << "number of time steps : " << timeStepNbr_ << "\n";
        ss << core::to_str(hybridModel_->state);
        return ss.str();
    }




    void initialize() override
    {
        if (integrator_ != nullptr)
            integrator_->initialize();
        else
            throw std::runtime_error("Error - Simulator has no integrator");
    }



    double startTime() override { return 0.; }
    double endTime() override { return 100.; }
    double timeStep() override { return 0.01; }

    void advance(double dt) override
    {
        if (integrator_)
            integrator_->advance(dt);
        else
            throw std::runtime_error("Error - no valid integrator in the simulator");
    }


    auto& getHybridModel() { return hybridModel_; }


    auto& getMHDModel() { return mhdModel_; }


    auto& getMultiPhysicsIntegrator() { return multiphysInteg_; }


private:
    auto find_model(std::string name)
    {
        return std::find(std::begin(modelNames_), std::end(modelNames_), name)
               != std::end(modelNames_);
    }

    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy_;
    std::unique_ptr<PHARE::amr::Integrator> integrator_;

    std::vector<std::string> modelNames_;
    std::vector<PHARE::amr::MessengerDescriptor> descriptors_;
    MessengerFactory messengerFactory_;

    float x_lo_[dimension];
    float x_up_[dimension];

    int maxLevelNumber_;
    double dt_;
    int timeStepNbr_;
    double finalTime_;

    // physical models that can be used
    std::shared_ptr<HybridModel> hybridModel_;
    std::shared_ptr<MHDModel> mhdModel_;

    std::shared_ptr<MultiPhysicsIntegrator> multiphysInteg_{nullptr};
};


struct SimulatorMaker
{
    SimulatorMaker(std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy)
        : hierarchy_{hierarchy}
    {
    }

    std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy_;

    template<typename Dimension, typename InterpOrder>
    std::unique_ptr<ISimulator> operator()(std::size_t userDim, std::size_t userInterpOrder,
                                           Dimension dimension, InterpOrder interp_order)
    {
        if (userDim == dimension() and userInterpOrder == interp_order())
        {
            std::size_t constexpr d  = dimension();
            std::size_t constexpr io = interp_order();

            PHARE::initializer::PHAREDict& theDict
                = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
            return std::make_unique<Simulator<d, io>>(theDict, hierarchy_);
        }
        else
        {
            return nullptr;
        }
    }
};


std::unique_ptr<PHARE::ISimulator> getSimulator(std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy)
{
    PHARE::initializer::PHAREDict& theDict
        = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    auto dim         = theDict["simulation"]["dimension"].template to<int>();
    auto interpOrder = theDict["simulation"]["interp_order"].template to<int>();
    return core::makeAtRuntime<SimulatorMaker>(dim, interpOrder, SimulatorMaker{hierarchy});
}



} // namespace PHARE
#endif
