#ifndef PHARE_SIMULATOR_SIMULATOR_HPP
#define PHARE_SIMULATOR_SIMULATOR_HPP



#include "phare_core.hpp"
#include "phare_types.hpp"

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/utilities/timestamps.hpp"
#include "amr/tagging/tagger_factory.hpp"

#include <chrono>
#include <exception>
#include <unordered_set>


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
    virtual bool dump(double timestamp, double timestep) { return false; } // overriding optional
};

template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
class Simulator : public ISimulator
{
public:
    NO_DISCARD double startTime() override { return startTime_; }
    NO_DISCARD double endTime() override { return finalTime_; }
    NO_DISCARD double timeStep() override { return dt_; }
    NO_DISCARD double currentTime() override { return currentTime_; }

    void initialize() override;
    double advance(double dt) override;

    std::vector<int> const& domainBox() const override { return hierarchy_->domainBox(); }
    std::vector<double> const& cellWidth() const override { return hierarchy_->cellWidth(); }
    std::size_t interporder() const override { return interp_order; }

    NO_DISCARD auto& getHybridModel() { return hybridModel_; }
    NO_DISCARD auto& getMHDModel() { return mhdModel_; }
    NO_DISCARD auto& getMultiPhysicsIntegrator() { return multiphysInteg_; }

    NO_DISCARD std::string to_str() override;

    bool dump(double timestamp, double timestep) override
    {
        if (rMan)
        {
            rMan->dump(timestamp, timestep);
        }

        if (dMan)
            return dMan->dump(timestamp, timestep);

        return false;
    }

    Simulator(PHARE::initializer::PHAREDict const& dict,
              std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy);
    ~Simulator()
    {
        if (coutbuf != nullptr)
            std::cout.rdbuf(coutbuf);
    }

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

    using SimFunctorParams = typename core::PHARE_Sim_Types::SimFunctorParams;
    using SimFunctors      = typename core::PHARE_Sim_Types::SimulationFunctors;

    using Integrator = PHARE::amr::Integrator<dimension>;


private:
    auto find_model(std::string name);

    auto static log_file_name()
    {
        // ".log" directory is not created here, but in python if PHARE_LOG != "NONE"
        if (auto log = core::get_env("PHARE_LOG"))
        {
            if (log == "RANK_FILES")
                return ".log/" + std::to_string(core::mpi::rank()) + ".out";


            if (log == "DATETIME_FILES")
            {
                auto date_time = core::mpi::date_time();
                auto rank      = std::to_string(core::mpi::rank());
                auto size      = std::to_string(core::mpi::size());
                return ".log/" + date_time + "_" + rank + "_of_" + size + ".out";
            }

            if (log != "NONE")
                throw std::runtime_error(
                    "PHARE_LOG invalid type, valid keys are RANK_FILES/DATETIME_FILES/NONE");
        }

        return std::string{""}; // unused
    }

    std::ofstream log_out{log_file_name()};
    std::streambuf* coutbuf = nullptr;
    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy_;
    std::unique_ptr<Integrator> integrator_;

    std::vector<std::string> modelNames_;
    std::vector<PHARE::amr::MessengerDescriptor> descriptors_;
    MessengerFactory messengerFactory_;

    float x_lo_[dimension];
    float x_up_[dimension];
    int maxLevelNumber_;
    double dt_;
    int timeStepNbr_           = 0;
    double startTime_          = 0;
    double finalTime_          = 0;
    double currentTime_        = 0;
    bool isInitialized         = false;
    std::size_t fineDumpLvlMax = 0;

    // physical models that can be used
    std::shared_ptr<HybridModel> hybridModel_;
    std::shared_ptr<MHDModel> mhdModel_;

    std::unique_ptr<PHARE::core::ITimeStamper> timeStamper;
    std::unique_ptr<PHARE::diagnostic::IDiagnosticsManager> dMan;
    std::unique_ptr<PHARE::restarts::IRestartsManager> rMan;

    SimFunctors functors_;

    SimFunctors functors_setup(PHARE::initializer::PHAREDict const& dict)
    {
        return {{"pre_advance", {/*empty vector*/}}};
    }

    std::shared_ptr<MultiPhysicsIntegrator> multiphysInteg_{nullptr};



    double restarts_init(initializer::PHAREDict const&);
    void diagnostics_init(initializer::PHAREDict const&);
    void hybrid_init(initializer::PHAREDict const&);
};



namespace
{
    inline auto logging(std::ofstream& log_out)
    {
        std::streambuf* buf = nullptr;
        if (auto log = core::get_env("PHARE_LOG"); log != "NONE")
        {
            buf = std::cout.rdbuf();
            std::cout.rdbuf(log_out.rdbuf());
        }
        return buf;
    }
} // namespace



//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------

template<std::size_t dim, std::size_t _interp, std::size_t nbRefinedPart>
double Simulator<dim, _interp, nbRefinedPart>::restarts_init(initializer::PHAREDict const& dict)
{
    rMan = restarts::RestartsManagerResolver::make_unique(*hierarchy_, *hybridModel_, dict);

    if (dict.contains("restart_time"))
        return (currentTime_ = dict["restart_time"].template to<double>());

    return 0;
}



template<std::size_t dim, std::size_t _interp, std::size_t nbRefinedPart>
void Simulator<dim, _interp, nbRefinedPart>::diagnostics_init(initializer::PHAREDict const& dict)
{
    dMan = PHARE::diagnostic::DiagnosticsManagerResolver::make_unique(*hierarchy_, *hybridModel_,
                                                                      dict);

    if (dict.contains("fine_dump_lvl_max"))
    {
        auto fine_dump_lvl_max = dict["fine_dump_lvl_max"].template to<int>();

        if (fine_dump_lvl_max > 0)
        { // copy for later
            this->fineDumpLvlMax                  = static_cast<std::size_t>(fine_dump_lvl_max);
            functors_["pre_advance"]["fine_dump"] = [&](SimFunctorParams const& params) {
                std::size_t level_nbr = params["level_nbr"].template to<int>();
                auto timestamp        = params["timestamp"].template to<double>();

                if (this->fineDumpLvlMax >= level_nbr)
                    this->dMan->dump_level(level_nbr, timestamp);
            };
        }
    }
}



template<std::size_t dim, std::size_t _interp, std::size_t nbRefinedPart>
void Simulator<dim, _interp, nbRefinedPart>::hybrid_init(initializer::PHAREDict const& dict)
{
    hybridModel_ = std::make_shared<HybridModel>(
        dict["simulation"], std::make_shared<typename HybridModel::resources_manager_type>());


    hybridModel_->resourcesManager->registerResources(hybridModel_->state);

    // we register the hybrid model for all possible levels in the hierarchy
    // since for now it is the only model available, same for the solver
    multiphysInteg_->registerModel(0, maxLevelNumber_ - 1, hybridModel_);

    multiphysInteg_->registerAndInitSolver(0, maxLevelNumber_ - 1,
                                           std::make_unique<SolverPPC>(dict["simulation"]["algo"]));

    multiphysInteg_->registerAndSetupMessengers(messengerFactory_);

    // hard coded for now, should get some params later from the dict
    auto hybridTagger_ = amr::TaggerFactory<PHARETypes>::make("HybridModel", "default");
    multiphysInteg_->registerTagger(0, maxLevelNumber_ - 1, std::move(hybridTagger_));

    if (dict["simulation"].contains("restarts"))
        startTime_ = restarts_init(dict["simulation"]["restarts"]);

    integrator_ = std::make_unique<Integrator>(dict, hierarchy_, multiphysInteg_, multiphysInteg_,
                                               startTime_, finalTime_);

    timeStamper = core::TimeStamperFactory::create(dict["simulation"]);

    if (dict["simulation"].contains("diagnostics"))
        diagnostics_init(dict["simulation"]["diagnostics"]);
}



template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
Simulator<_dimension, _interp_order, _nbRefinedPart>::Simulator(
    PHARE::initializer::PHAREDict const& dict,
    std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy)
    : coutbuf{logging(log_out)}
    , hierarchy_{hierarchy}
    , modelNames_{"HybridModel"}
    , descriptors_{PHARE::amr::makeDescriptors(modelNames_)}
    , messengerFactory_{descriptors_}
    , maxLevelNumber_{dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>()}
    , dt_{dict["simulation"]["time_step"].template to<double>()}
    , timeStepNbr_{dict["simulation"]["time_step_nbr"].template to<int>()}
    , finalTime_{dt_ * timeStepNbr_}
    , functors_{functors_setup(dict)}
    , multiphysInteg_{std::make_shared<MultiPhysicsIntegrator>(dict["simulation"], functors_)}
{
    if (find_model("HybridModel"))
        hybrid_init(dict);
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

    if (core::mpi::any(core::Errors::instance().any()))
    {
        this->dMan.release(); // closes/flushes hdf5 files
        throw std::runtime_error("forcing error");
    }

    isInitialized = true;

    if (hierarchy_->isFromRestart())
        hierarchy_->closeRestartFile();
}




template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
double Simulator<_dimension, _interp_order, _nbRefinedPart>::advance(double dt)
{
    double dt_new = 0;

    if (!integrator_)
        throw std::runtime_error("Error - no valid integrator in the simulator");

    try
    {
        PHARE_LOG_SCOPE(1, "Simulator::advance");
        dt_new       = integrator_->advance(dt);
        currentTime_ = startTime_ + ((*timeStamper) += dt);
    }
    catch (std::runtime_error const& e)
    {
        std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
        std::rethrow_exception(std::current_exception());
    }
    catch (...)
    {
        std::cerr << "UNKNOWN EXCEPTION CAUGHT" << std::endl;
        std::rethrow_exception(std::current_exception());
    }

    if (core::mpi::any(core::Errors::instance().any()))
    {
        this->dMan.release(); // closes/flushes hdf5 files
        throw std::runtime_error("forcing error");
    }

    return dt_new;
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


} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_H*/
