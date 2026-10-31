#ifndef PHARE_SIMULATOR_SIMULATOR_DECL_HPP
#define PHARE_SIMULATOR_SIMULATOR_DECL_HPP

// Simulator class and its light members only. The heavy members (those that instantiate a
// concrete model/solver/messenger/...) are defined in the simulator_<part>.hpp headers, so
// each part can be explicitly instantiated in its own translation unit to bound compiler
// memory, see src/python3/cpp_simulator_*.cpp. simulator/simulator.hpp includes every part.

#include "phare_solver.hpp"

#include "core/def.hpp"
#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/thread_pool.hpp"
#include "mpi/mpi_utils.hpp"
#include "core/utilities/timestamps.hpp"

#include "amr/wrappers/hierarchy.hpp"
#include "amr/wrappers/integrator.hpp"
#include "amr/messengers/refinement_config.hpp"
#include "amr/messengers/messenger_factory.hpp"

#include "diagnostic/diagnostic_manager.hpp"

#include "restarts/restarts_manager.hpp"
#include "simulator/simulator_def.hpp"

#include <stdexcept>
#include <vector>
#include <string>
#include <optional>
#include <algorithm>


namespace PHARE
{

template<auto opts>
class Simulator : public ISimulator
{
public:
    auto static constexpr options = opts;

    std::size_t static constexpr dimension     = opts.dimension;
    std::size_t static constexpr interp_order  = opts.interp_order;
    std::size_t static constexpr nbRefinedPart = opts.nbRefinedPart;

    using SAMRAITypes            = PHARE::amr::SAMRAI_Types;
    using PHARETypes             = solver::PHARE_Types<opts>;
    using IPhysicalModel         = PHARE::solver::IPhysicalModel<SAMRAITypes>;
    using MessengerFactory       = PHARETypes::MessengerFactory;
    using MultiPhysicsIntegrator = PHARETypes::MultiPhysicsIntegrator_t;
    using SimFunctorParams       = core::PHARE_Sim_Types::SimFunctorParams;
    using SimFunctors            = core::PHARE_Sim_Types::SimulationFunctors;
    using Integrator             = PHARE::amr::Integrator<dimension>;


    // does nothing that can throw, see setup()
    Simulator(std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy) noexcept
        : hierarchy_{hierarchy}
    {
    }

    // throwing setup, the destructor still runs if it throws
    void setup(PHARE::initializer::PHAREDict const& dict);

    // the thread pools are configured in setup(), for this simulation only
    ~Simulator() { core::ThreadPool::INSTANCE().stop(); }


    NO_DISCARD double startTime() override { return startTime_; }
    NO_DISCARD double endTime() override { return finalTime_; }
    NO_DISCARD double timeStep() override { return dt_; }
    NO_DISCARD double currentTime() override { return currentTime_; }

    void initialize() override;
    double advance() override;

    std::vector<int> const& domainBox() const override { return hierarchy_->domainBox(); }
    std::vector<double> const& cellWidth() const override { return hierarchy_->cellWidth(); }
    std::size_t interporder() const override { return interp_order; }

    NO_DISCARD auto& getHybridModel() { return hyb_.model_; }
    NO_DISCARD auto& getMHDModel() { return mhd_.model_; }
    NO_DISCARD auto& getMultiPhysicsIntegrator() { return multiphysInteg_; }

    NO_DISCARD std::string to_str() const override;

    bool dump_diagnostics(double timestamp, double timestep) override
    {
        if (dMan)
            return dMan->dump(timestamp, timestep);
        return false;
    }
    bool dump_restarts(double timestamp, double timestep) override
    {
        if (rMan)
            return rMan->dump(timestamp, timestep);
        return false;
    }

    std::string summary() const;




protected:
    // provided to force flush for diags
    void reset_dman() { this->dMan.reset(); }

private:
    auto find_model(std::string name);

    std::unique_ptr<std::ofstream> static log_file();

    std::optional<CoutRedirect> coutRedirect_;
    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy_;
    std::unique_ptr<Integrator> integrator_;

    std::vector<std::string> modelNames_;
    std::vector<PHARE::amr::MessengerDescriptor> descriptors_;
    MessengerFactory messengerFactory_{{}};


    int maxLevelNumber_              = 0;
    int maxMHDLevel_                 = 0;
    double dt_                       = 0;
    core::TimeStepType timeStepType_ = core::TimeStepType::constant;
    double cflWave_                  = 0;
    double cflDiffusive_             = 0;
    int timeStepNbr_                 = 0;
    double startTime_                = 0;
    double finalTime_                = 0;
    double currentTime_              = 0;
    std::size_t fineDumpLvlMax       = 0;
    bool isInitialized               = false;
    bool allowEmergencyDumps         = false;

    // empty when opts has the corresponding model disabled, see phare_solver.hpp
    solver::HybridSimState<opts> hyb_;
    solver::MHDSimState<opts> mhd_;

    std::unique_ptr<PHARE::core::ITimeStamper> timeStamper;
    std::unique_ptr<PHARE::diagnostic::IDiagnosticsManager> dMan;
    std::unique_ptr<PHARE::restarts::IRestartsManager> rMan;

    SimFunctors functors_;

    SimFunctors functors_setup(PHARE::initializer::PHAREDict const& dict)
    {
        return {{"pre_advance", {/*empty vector*/}}};
    }

    std::shared_ptr<MultiPhysicsIntegrator> multiphysInteg_{nullptr};



    double restart_time(initializer::PHAREDict const&);
    void threading_init(initializer::PHAREDict const&);
    void hybrid_init(initializer::PHAREDict const&);
    void mhd_init(initializer::PHAREDict const&);

    // per part definitions, bodies are empty if opts disables the model

    // simulator_models.hpp
    void hybrid_model_init(initializer::PHAREDict const&);
    void mhd_model_init(initializer::PHAREDict const&);
    std::string models_to_str() const;

    // simulator_integrator.hpp
    void multiphysics_init(initializer::PHAREDict const&);
    void hybrid_register_model();
    void mhd_register_model();
    void hybrid_tagger_init(initializer::PHAREDict const&);
    void mhd_tagger_init(initializer::PHAREDict const&);
    void finalize_init(PHARE::initializer::PHAREDict const& dict);

    // simulator_solvers.hpp
    void hybrid_solver_init(initializer::PHAREDict const&);
    void mhd_solver_init(initializer::PHAREDict const&);

    // simulator_messengers.hpp
    void messengers_init();

    // simulator_diagnostics.hpp
    void diagnostics_init(initializer::PHAREDict const&, auto&);
    void hybrid_diagnostics_init(initializer::PHAREDict const&);
    void mhd_diagnostics_init(initializer::PHAREDict const&);

    // simulator_restarts.hpp
    void hybrid_restarts_init(initializer::PHAREDict const&);
    void mhd_restarts_init(initializer::PHAREDict const&);

    void handle_dictionary_exception(core::DictionaryException const& ex);
};

template<auto opts>
std::unique_ptr<std::ofstream> Simulator<opts>::log_file()
{
    using namespace PHARE::core;

    std::string const base = ".log/";
    std::string const ext  = ".out";

    // ".log" directory is not created here, but in simulator.py
    if (auto log = get_env("PHARE_LOG"))
    {
        if (log == "RANK_FILES")
            return std::make_unique<std::ofstream>(base + std::to_string(mpi::rank()) + ext);

        if (log == "DATETIME_FILES")
        {
            auto date_time = mpi::date_time();
            auto rank      = std::to_string(mpi::rank());
            auto size      = std::to_string(mpi::size());
            return std::make_unique<std::ofstream>(base + date_time + "_" + rank + "_of_" + size
                                                   + ext);
        }

        if (log == "NULL")
            return std::make_unique<std::ofstream>("/dev/null");

        if (log != "CLI")
            throw std::runtime_error(
                "PHARE_LOG invalid type, valid keys are RANK_FILES/DATETIME_FILES/CLI/NULL");
    }

    return nullptr;
}

//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------

template<auto opts>
double Simulator<opts>::restart_time(initializer::PHAREDict const& dict)
{
    return cppdict::get_value(dict, "simulation/restarts/restart_time", 0.);
}



// env vars override the input, e.g. for benchmarks (tools/bench/functional/permutor.py)
template<auto opts>
void Simulator<opts>::threading_init(initializer::PHAREDict const& dict)
{
    auto const input = [&](std::string const& key) {
        return static_cast<std::size_t>(cppdict::get_value(dict, "simulation/threading/" + key, 1));
    };
    auto const pools   = core::get_env_as("PHARE_THREAD_POOLS", input("pools"));
    auto const threads = core::get_env_as("PHARE_THREADS_PER_POOL", input("threads_per_pool"));
    core::ThreadPool::INSTANCE().configure(pools, threads);
}


template<auto opts>
void Simulator<opts>::hybrid_init(initializer::PHAREDict const& dict)
{
    hybrid_model_init(dict);

    // we register the hybrid model for all possible levels in the hierarchy
    // since for now it is the only model available, same for the solver
    hybrid_register_model();
    hybrid_solver_init(dict);
    messengers_init();
    hybrid_tagger_init(dict);

    if (dict["simulation"].contains("diagnostics"))
        hybrid_diagnostics_init(dict["simulation"]["diagnostics"]);
}


template<auto opts>
void Simulator<opts>::mhd_init(initializer::PHAREDict const& dict)
{
    mhd_model_init(dict);

    // we register the mhd model for all possible levels in the hierarchy
    // since for now it is the only model available, same for the solver
    mhd_register_model();
    mhd_solver_init(dict);
    messengers_init();
    mhd_tagger_init(dict);

    if (dict["simulation"].contains("diagnostics"))
        mhd_diagnostics_init(dict["simulation"]["diagnostics"]);
}



template<auto opts>
void Simulator<opts>::setup(PHARE::initializer::PHAREDict const& dict)
{
    if (bool const already_setup = multiphysInteg_ != nullptr; already_setup)
        throw std::runtime_error("Simulator::setup already called");

    coutRedirect_.emplace(log_file());
    threading_init(dict);

    if (!hierarchy_)
        throw std::runtime_error("NO HIERARCHY!");

    modelNames_       = dict["simulation"]["models"].template to<std::vector<std::string>>();
    descriptors_      = PHARE::amr::makeDescriptors(modelNames_);
    messengerFactory_ = MessengerFactory{descriptors_, PHARE::amr::RefinementConfig::FROM(dict)};
    maxLevelNumber_   = dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>();
    maxMHDLevel_      = dict["simulation"]["AMR"]["max_mhd_level"].template to<int>();
    dt_               = cppdict::get_value(dict, "simulation/time_step/value", 0.);
    timeStepType_
        = cppdict::get_value(dict, "simulation/time_step/mode", core::TimeStepType::constant);
    cflWave_      = cppdict::get_value(dict, "simulation/time_step/cfl_wave", 0.);
    cflDiffusive_ = cppdict::get_value(dict, "simulation/time_step/cfl_diffusive", 0.);
    timeStepNbr_  = cppdict::get_value(dict, "simulation/time_step_nbr", 0);
    finalTime_    = dict["simulation"]["final_time"].template to<double>();
    functors_     = functors_setup(dict);
    multiphysics_init(dict);

    currentTime_ = restart_time(dict);
    // finalTime_ is computed in Python as start_time + time_step_nbr * time_step
    // (start_time == restart_time), so it already accounts for the restart offset.

    // we would need a different restart manager for mhd and hybrid if both models are used

    // the build enabling a model (static_assert in solver::PHARE_Types) is not enough: the dict
    // still has to ask for one of them by name.
    bool anyModel = false;

    if constexpr (has_hybrid_v<opts>)
    {
        if (find_model("HybridModel"))
        {
            anyModel = true;
            hybrid_init(dict);
            if (dict["simulation"].contains("restarts"))
                hybrid_restarts_init(dict["simulation"]["restarts"]);
        }
    }
    else
    {
        if (find_model("HybridModel"))
            throw std::runtime_error("HybridModel requested but this build has no hybrid support");
    }

    if constexpr (has_mhd_v<opts>)
    {
        if (find_model("MHDModel"))
        {
            anyModel = true;
            mhd_init(dict);
            if (dict["simulation"].contains("restarts"))
                mhd_restarts_init(dict["simulation"]["restarts"]);
        }
    }
    else
    {
        if (find_model("MHDModel"))
            throw std::runtime_error("MHDModel requested but this build has no MHD support");
    }

    if (!anyModel)
    {
        std::string names;
        for (auto const& name : modelNames_)
            names += (names.empty() ? "" : ", ") + name;
        throw std::runtime_error("unsupported model, none of [" + names
                                 + "] is supported by this build");
    }

    finalize_init(dict);

    amr::ResourcesManagerGlobals::registerForRestarts();
}



template<auto opts>
std::string Simulator<opts>::to_str() const
{
    std::stringstream ss;
    ss << "PHARE SIMULATOR\n";
    ss << "------------------------------------\n";
    ss << "interpolation order  : " << interp_order << "\n";
    ss << "dimension            : " << dimension << "\n";
    ss << "time step            : " << dt_ << "\n";
    ss << "number of time steps : " << timeStepNbr_ << "\n";
    ss << models_to_str();
    return ss.str();
}

template<auto opts>
void Simulator<opts>::initialize()
{
    PHARE_LOG_SCOPE(1, "Simulator::initialize");

    if (bool const needs_setup = multiphysInteg_ == nullptr; needs_setup)
        throw std::runtime_error("Simulator::initialize: setup() needs to be called first");

    std::optional<std::string> error = std::nullopt;

    try
    {
        if (isInitialized)
            throw std::runtime_error("cannot initialize  - simulator already isInitialized");

        if (integrator_ == nullptr)
            throw std::runtime_error("Error - Simulator has no integrator");

        integrator_->initialize();

        // First dt_ for the init dump
        if (timeStepType_ == core::TimeStepType::adaptive)
        {
            double const dt = multiphysInteg_->computeStableDt(
                *hierarchy_, solver::CFLNumbers{cflWave_, cflDiffusive_});
            dt_ = std::min(dt, finalTime_ - currentTime_);
        }
    }
    catch (core::DictionaryException const& ex)
    {
        handle_dictionary_exception(ex);
        error = std::string{"DICTIONARY EXCEPTION CAUGHT: "} + ex.what();
        PHARE_LOG_ERROR(*error);
    }
    catch (std::exception const& e)
    {
        error = std::string{"EXCEPTION CAUGHT: "} + e.what();
        PHARE_LOG_ERROR(*error);
    }
    catch (...)
    {
        error = "UNKNOWN EXCEPTION CAUGHT";
        PHARE_LOG_ERROR(*error);
    }

    if (mpi::any_errors())
    {
        this->dMan.reset(); // closes/flushes hdf5 files
        if (error)
            throw std::runtime_error(*error);
        throw std::runtime_error("forcing error");
    }

    isInitialized = true;

    if (hierarchy_->isFromRestart())
        hierarchy_->closeRestartFile();
}


template<auto opts>
double Simulator<opts>::advance()
{
    PHARE_LOG_SCOPE(1, "Simulator::advance");

    std::optional<std::string> error = std::nullopt;

    try
    {
        if (!integrator_)
            throw std::runtime_error("Error - no valid integrator in the simulator");

        if (timeStepType_ == core::TimeStepType::adaptive)
        {
            double const dt = multiphysInteg_->computeStableDt(
                *hierarchy_, solver::CFLNumbers{cflWave_, cflDiffusive_});
            dt_ = std::min(dt, finalTime_ - currentTime_);
        }

        integrator_->advance(dt_);
        currentTime_ = startTime_ + ((*timeStamper) += dt_);
    }
    catch (core::DictionaryException const& ex)
    {
        handle_dictionary_exception(ex);
        error = std::string{"DICTIONARY EXCEPTION CAUGHT: "} + ex.what();
        PHARE_LOG_ERROR(*error);
    }
    catch (std::exception const& e)
    {
        error = std::string{"EXCEPTION CAUGHT: "} + e.what();
        PHARE_LOG_ERROR(*error);
    }
    catch (...)
    {
        error = "UNKNOWN EXCEPTION CAUGHT";
        PHARE_LOG_ERROR(*error);
    }

    if (mpi::any_errors())
    {
        this->dMan.reset(); // closes/flushes hdf5 files
        if (error)
            throw std::runtime_error(*error);
        throw std::runtime_error("forcing error");
    }

    return dt_;
}

template<auto opts>
void Simulator<opts>::handle_dictionary_exception(core::DictionaryException const& ex)
{
    constexpr static std::array dump_exceptions{"Updater::updatePopulations",
                                                "HybridLevelInitializer::initialize",
                                                "SolverMHD::advanceLevel"};
    if (this->allowEmergencyDumps)
        for (auto const& exception_id : dump_exceptions)
            if (ex.id() == exception_id)
                for (int ilvl = 0; ilvl < hierarchy_->getMaxNumberOfLevels(); ++ilvl)
                    this->dMan->dump_level(ilvl, currentTime_);
}

template<auto opts>
auto Simulator<opts>::find_model(std::string name)
{
    if (modelNames_.empty())
        throw std::runtime_error("Simulator: No models found!");
    return std::find(std::begin(modelNames_), std::end(modelNames_), name) != std::end(modelNames_);
}


template<typename Simulator>
std::unique_ptr<Simulator> makeSimulator(std::shared_ptr<amr::Hierarchy> const& hierarchy)
{
    auto sim = std::make_unique<Simulator>(hierarchy);
    sim->setup(initializer::PHAREDictHandler::INSTANCE().dict()); // owned: destructed on throw
    return sim;
}



} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_DECL_HPP*/
