
#ifndef PHARE_SIMULATOR_SIMULATOR_H
#define PHARE_SIMULATOR_SIMULATOR_H

#include "include.h"
#include "phare_types.h"


namespace PHARE
{
class ISimulator
{
public:
    virtual void initialize()               = 0;
    virtual double startTime()              = 0;
    virtual double endTime()                = 0;
    virtual double currentTime()            = 0;
    virtual double timeStep()               = 0;
    virtual void advance(double dt)         = 0;
    virtual std::string to_str()            = 0;
    virtual std::string domainBox() const   = 0;
    virtual std::string cellWidth() const   = 0;
    virtual std::size_t interporder() const = 0;
    virtual ~ISimulator() {}
};



template<std::size_t _dimension, std::size_t _interp_order, size_t _nbRefinedPart>
class Simulator : public ISimulator
{
public:
    static constexpr size_t dimension     = _dimension;
    static constexpr size_t interp_order  = _interp_order;
    static constexpr size_t nbRefinedPart = _nbRefinedPart;

    using SAMRAITypes = PHARE::amr::SAMRAI_Types;
    using PHARETypes  = PHARE_Types<dimension, interp_order, nbRefinedPart>;

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


            integrator_ = std::make_unique<PHARE::amr::DimIntegrator<dimension>>(
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



    double startTime() override { return 0.; }

    double endTime() override { return finalTime_; }
    double timeStep() override { return dt_; }
    double currentTime() override { return currentTime_; }

    std::string domainBox() const override { return hierarchy_->domainBox(); }
    std::string cellWidth() const override { return hierarchy_->cellWidth(); }
    std::size_t interporder() const override { return interp_order; }

    void advance(double dt) override
    {
        try
        {
            if (integrator_)
            {
                integrator_->advance(dt);
                currentTime_ += dt;
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
    int timeStepNbr_    = 0;
    double finalTime_   = 0;
    double currentTime_ = 0;
    bool isInitialized  = false;

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

    template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
    std::unique_ptr<ISimulator> operator()(std::size_t userDim, std::size_t userInterpOrder,
                                           size_t userNbRefinedPart, Dimension dimension,
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


} // namespace PHARE
#endif
