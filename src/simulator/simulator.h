
#ifndef PHARE_SIMULATOR_SIMULATOR_H
#define PHARE_SIMULATOR_SIMULATOR_H

#include "include.h"
#include "phare_types.h"
#include <chrono>

namespace PHARE
{
class ISimulator
{
public:
    virtual double startTime()   = 0;
    virtual double endTime()     = 0;
    virtual double currentTime() = 0;
    virtual double timeStep()    = 0;

    virtual void initialize()         = 0;
    virtual double advance(double dt) = 0;

    virtual std::vector<int> const& domainBox() const    = 0;
    virtual std::vector<double> const& cellWidth() const = 0;
    virtual std::size_t interporder() const              = 0;

    virtual std::string to_str() = 0;

    virtual ~ISimulator() {}
};



template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
class Simulator : public ISimulator
{
public:
    double startTime() override { return 0.; }
    double endTime() override { return finalTime_; }
    double timeStep() override { return dt_; }
    double currentTime() override { return currentTime_; }

    void initialize() override;
    double advance(double dt) override;

    std::vector<int> const& domainBox() const override { return hierarchy_->domainBox(); }
    std::vector<double> const& cellWidth() const override { return hierarchy_->cellWidth(); }
    std::size_t interporder() const override { return interp_order; }

    auto& getHybridModel() { return hybridModel_; }
    auto& getMHDModel() { return mhdModel_; }
    auto& getMultiPhysicsIntegrator() { return multiphysInteg_; }

    std::string to_str() override;

    Simulator(PHARE::initializer::PHAREDict dict,
              std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy);


    static constexpr std::size_t dimension     = _dimension;
    static constexpr std::size_t interp_order  = _interp_order;
    static constexpr std::size_t nbRefinedPart = _nbRefinedPart;

    using SAMRAITypes = PHARE::amr::SAMRAI_Types;
    using PHARETypes  = PHARE_Types<dimension, interp_order, nbRefinedPart>;

    using IPhysicalModel = PHARE::solver::IPhysicalModel<SAMRAITypes>;
    using HybridModel    = typename PHARETypes::HybridModel_t;
    using MHDModel       = typename PHARETypes::MHDModel_t;

    using SolverMHD = typename PHARETypes::SolverMHD_t;
    using SolverPPC = typename PHARETypes::SolverPPC_t;

    using MessengerFactory       = typename PHARETypes::MessengerFactory;
    using MultiPhysicsIntegrator = typename PHARETypes::MultiPhysicsIntegrator;

    using Integrator = PHARE::amr::Integrator<dimension>;


private:
    auto find_model(std::string name);

    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy_;
    std::unique_ptr<Integrator> integrator_;

    std::vector<std::string> modelNames_;
    std::vector<PHARE::amr::MessengerDescriptor> descriptors_;
    MessengerFactory messengerFactory_;

    float x_lo_[dimension];
    float x_up_[dimension];
    int maxLevelNumber_;
    double dt_;
    int timeStepNbr_    = 0;
    double finalTime_   = 0;
    double currentTime_ = 0;
    bool isInitialized  = false;

    // physical models that can be used
    std::shared_ptr<HybridModel> hybridModel_;
    std::shared_ptr<MHDModel> mhdModel_;

    std::shared_ptr<MultiPhysicsIntegrator> multiphysInteg_{nullptr};
};


//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------



template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
Simulator<_dimension, _interp_order, _nbRefinedPart>::Simulator(
    PHARE::initializer::PHAREDict dict, std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy)
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
            dict["simulation"], std::make_shared<typename HybridModel::resources_manager_type>());


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


        integrator_ = std::make_unique<Integrator>(dict, hierarchy, multiphysInteg_,
                                                   multiphysInteg_, startTime, endTime);
    }
    else
        throw std::runtime_error("unsupported model");
}



template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
std::string Simulator<_dimension, _interp_order, _nbRefinedPart>::to_str()
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




template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
void Simulator<_dimension, _interp_order, _nbRefinedPart>::initialize()
{
    try
    {
        if (isInitialized)
            std::runtime_error("cannot initialize  - simulator already isInitialized");

        if (integrator_ != nullptr)
            integrator_->initialize();
        else
            throw std::runtime_error("Error - Simulator has no integrator");
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
        std::rethrow_exception(std::current_exception());
    }
    catch (...)
    {
        std::cerr << "UNKNOWN EXCEPTION CAUGHT" << std::endl;
        std::rethrow_exception(std::current_exception());
    }
    isInitialized = true;
}




template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
double Simulator<_dimension, _interp_order, _nbRefinedPart>::advance(double dt)
{
    try
    {
        if (integrator_)
        {
            auto dt_new = integrator_->advance(dt);
            currentTime_ += dt;
            return dt_new;
        }
        else
            throw std::runtime_error("Error - no valid integrator in the simulator");
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
        std::rethrow_exception(std::current_exception());
    }
    catch (...)
    {
        std::cerr << "UNKNOWN EXCEPTION CAUGHT" << std::endl;
        std::rethrow_exception(std::current_exception());
    }
}




template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
auto Simulator<_dimension, _interp_order, _nbRefinedPart>::find_model(std::string name)
{
    return std::find(std::begin(modelNames_), std::end(modelNames_), name) != std::end(modelNames_);
}



struct SimulatorMaker
{
    SimulatorMaker(std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy)
        : hierarchy_{hierarchy}
    {
    }

    std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy_;

    template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
    std::unique_ptr<ISimulator> operator()(std::size_t userDim, std::size_t userInterpOrder,
                                           std::size_t userNbRefinedPart, Dimension dimension,
                                           InterpOrder interp_order, NbRefinedPart nbRefinedPart)
    {
        if (userDim == dimension() and userInterpOrder == interp_order()
            and userNbRefinedPart == nbRefinedPart())
        {
            std::size_t constexpr d  = dimension();
            std::size_t constexpr io = interp_order();
            std::size_t constexpr nb = nbRefinedPart();

            PHARE::initializer::PHAREDict& theDict
                = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
            return std::make_unique<Simulator<d, io, nb>>(theDict, hierarchy_);
        }
        else
        {
            return nullptr;
        }
    }
};


std::unique_ptr<PHARE::ISimulator> getSimulator(std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy);


template<std::size_t dim, std::size_t interp, size_t nbRefinedPart>
std::unique_ptr<Simulator<dim, interp, nbRefinedPart>>
makeSimulator(std::shared_ptr<amr::Hierarchy> const& hierarchy)
{
    return std::make_unique<Simulator<dim, interp, nbRefinedPart>>(
        initializer::PHAREDictHandler::INSTANCE().dict(), hierarchy);
}



struct SimulatorDiagnostics
{
    SimulatorDiagnostics(PHARE::ISimulator& simulator, PHARE::amr::Hierarchy& hierarchy)
    {
        auto dict  = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
        this->dMan = PHARE::core::makeAtRuntime<Maker>(
            dict["simulation"]["dimension"].template to<int>(),
            dict["simulation"]["interp_order"].template to<int>(),
            dict["simulation"]["refined_particle_nbr"].template to<int>(),
            Maker{simulator, hierarchy});
        if (!this->dMan)
            throw std::runtime_error("Runtime diagnostic deduction failed");
    }

    struct Maker
    {
        Maker(PHARE::ISimulator& _simulator, PHARE::amr::Hierarchy& _hierarchy)
            : hierarchy{_hierarchy}
            , simulator{_simulator}
        {
        }


        template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
        std::unique_ptr<PHARE::diagnostic::IDiagnosticsManager>
        operator()(std::size_t userDim, std::size_t userInterpOrder, std::size_t userNbRefinedPart,
                   Dimension dimension, InterpOrder interp_order, NbRefinedPart nbRefinedPart)
        {
            auto& dict = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
            if (dict["simulation"].contains("diagnostics"))
            {
                if (userDim == dimension() and userInterpOrder == interp_order()
                    and userNbRefinedPart == nbRefinedPart())
                {
                    constexpr std::size_t d  = dimension();
                    constexpr std::size_t io = interp_order();
                    constexpr std::size_t nb = nbRefinedPart();

                    auto& cast_simulator = dynamic_cast<PHARE::Simulator<d, io, nb>&>(simulator);

                    return PHARE::diagnostic::DiagnosticsManagerResolver::make_unique(
                        hierarchy, *cast_simulator.getHybridModel(),
                        dict["simulation"]["diagnostics"]);
                }
            }
            return nullptr;
        }

        PHARE::amr::Hierarchy& hierarchy;
        PHARE::ISimulator& simulator;
    };

    void dump(double timestamp, double timestep) { dMan->dump(timestamp, timestep); }

    std::unique_ptr<PHARE::diagnostic::IDiagnosticsManager> dMan;
};


} // namespace PHARE
#endif
