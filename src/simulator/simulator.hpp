#ifndef PHARE_SIMULATOR_SIMULATOR_HPP
#define PHARE_SIMULATOR_SIMULATOR_HPP


#include "phare_core.hpp"
#include "phare_types.hpp"

#include "core/def.hpp"
#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/utilities/timestamps.hpp"

#include "amr/wrappers/integrator.hpp"
#include "amr/tagging/tagger_factory.hpp"
#include "amr/load_balancing/load_balancer_details.hpp"
#include "amr/load_balancing/load_balancer_manager.hpp"
#include "amr/load_balancing/load_balancer_estimator_hybrid.hpp"
#include "amr/load_balancing/load_balancer_estimator_mhd.hpp"

#include "python3/mhd_defaults/default_mhd_time_stepper.hpp"

#include "diagnostic/diagnostics.hpp"

#include "restarts/restarts.hpp"

#include <stdexcept>
#include <vector>
#include <string>


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

template<auto opts>
class Simulator : public ISimulator
{
public:
    std::size_t static constexpr dimension     = opts.dimension;
    std::size_t static constexpr interp_order  = opts.interp_order;
    std::size_t static constexpr nbRefinedPart = opts.nbRefinedPart;

    using SAMRAITypes            = PHARE::amr::SAMRAI_Types;
    using PHARETypes             = PHARE_Types<opts>;
    using IPhysicalModel         = PHARE::solver::IPhysicalModel<SAMRAITypes>;
    using HybridModel            = PHARETypes::HybridModel_t;
    using MHDModel               = PHARETypes::MHDModel_t;
    using SolverMHD              = PHARETypes::SolverMHD_t;
    using SolverPPC              = PHARETypes::SolverPPC_t;
    using MessengerFactory       = PHARETypes::MessengerFactory;
    using MultiPhysicsIntegrator = PHARETypes::MultiPhysicsIntegrator;
    using SimFunctorParams       = core::PHARE_Sim_Types::SimFunctorParams;
    using SimFunctors            = core::PHARE_Sim_Types::SimulationFunctors;
    using Integrator             = PHARE::amr::Integrator<dimension>;

    using HybridResourceManager_t = HybridModel::resources_manager_type;
    using MHDResourceManager_t    = MHDModel::resources_manager_type;


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
        {
            return dMan->dump(timestamp, timestep);
        }

        return false;
    }

    Simulator(PHARE::initializer::PHAREDict const& dict,
              std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy);
    ~Simulator()
    {
        if (coutbuf != nullptr)
            std::cout.rdbuf(coutbuf);
    }



protected:
    // provided to force flush for diags
    void reset_dman() { this->dMan.reset(); }

private:
    auto find_model(std::string name);

    std::unique_ptr<std::ofstream> static log_file()
    {
        // ".log" directory is not created here, but in simulator.py
        if (auto log = core::get_env("PHARE_LOG"))
        {
            if (log == "RANK_FILES")
                return std::make_unique<std::ofstream>(".log/" + std::to_string(core::mpi::rank())
                                                       + ".out");


            if (log == "DATETIME_FILES")
            {
                auto date_time = core::mpi::date_time();
                auto rank      = std::to_string(core::mpi::rank());
                auto size      = std::to_string(core::mpi::size());
                return std::make_unique<std::ofstream>(".log/" + date_time + "_" + rank + "_of_"
                                                       + size + ".out");
            }

            if (log == "NULL")
                return std::make_unique<std::ofstream>("/dev/null");

            if (log != "CLI")
                throw std::runtime_error(
                    "PHARE_LOG invalid type, valid keys are RANK_FILES/DATETIME_FILES/CLI/NULL");
        }

        return nullptr;
    }

    std::unique_ptr<std::ofstream> log_out{log_file()};
    std::streambuf* coutbuf = nullptr;
    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy_;
    std::unique_ptr<Integrator> integrator_;

    std::vector<std::string> modelNames_;
    std::vector<PHARE::amr::MessengerDescriptor> descriptors_;
    MessengerFactory messengerFactory_;

    float x_lo_[dimension];
    float x_up_[dimension];
    int maxLevelNumber_;
    int maxMHDLevel_;
    double dt_;
    int timeStepNbr_           = 0;
    double startTime_          = 0;
    double finalTime_          = 0;
    double currentTime_        = 0;
    bool isInitialized         = false;
    std::size_t fineDumpLvlMax = 0;

    std::shared_ptr<HybridResourceManager_t> hyb_resman_ptr;
    std::shared_ptr<MHDResourceManager_t> mhd_resman_ptr;

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



    double restart_time(initializer::PHAREDict const&);
    void diagnostics_init(initializer::PHAREDict const&, auto&);
    void hybrid_init(initializer::PHAREDict const&);
    void mhd_init(initializer::PHAREDict const&);
};



namespace
{
    inline auto logging(std::unique_ptr<std::ofstream>& log_out)
    {
        std::streambuf* buf = nullptr;
        if (log_out)
        {
            buf = std::cout.rdbuf();
            std::cout.rdbuf(log_out->rdbuf());
        }
        return buf;
    }
} // namespace



//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------

template<auto opts>
double Simulator<opts>::restart_time(initializer::PHAREDict const& dict)
{
    return cppdict::get_value(dict, "simulation/restarts/restart_time", 0.);
}



template<auto opts>
void Simulator<opts>::diagnostics_init(initializer::PHAREDict const& dict, auto& model)
{
    dMan = PHARE::diagnostic::DiagnosticsManagerResolver::make_unique(*hierarchy_, model, dict);

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


template<auto opts>
void Simulator<opts>::hybrid_init(initializer::PHAREDict const& dict)
{
    hybridModel_ = std::make_shared<HybridModel>(dict["simulation"], hyb_resman_ptr);
    hyb_resman_ptr->registerResources(hybridModel_->state); // still valid, never moved

    // we register the hybrid model for all possible levels in the hierarchy
    // since for now it is the only model available, same for the solver
    multiphysInteg_->registerModel(maxMHDLevel_, maxLevelNumber_ - 1, hybridModel_);

    multiphysInteg_->registerAndInitSolver(maxMHDLevel_, maxLevelNumber_ - 1,
                                           std::make_unique<SolverPPC>(dict["simulation"]["algo"]));

    multiphysInteg_->registerAndSetupMessengers(messengerFactory_);

    // hard coded for now, should get some params later from the dict
    if (dict["simulation"]["AMR"]["refinement"].contains("tagging"))
    {
        if (dict["simulation"]["AMR"]["refinement"]["tagging"]["method"].template to<std::string>()
            != "none")
        {
            auto hybridTagger_ = amr::TaggerFactory<HybridModel>::make(
                dict["simulation"]["AMR"]["refinement"]["tagging"]);
            multiphysInteg_->registerTagger(maxMHDLevel_, maxLevelNumber_ - 1,
                                            std::move(hybridTagger_));
        }
    }

    amr::LoadBalancerDetails lb_info
        = amr::LoadBalancerDetails::FROM(dict["simulation"]["AMR"]["loadbalancing"]);

    auto lbm_ = std::make_unique<amr::LoadBalancerManager<opts.dimension>>(dict);
    auto lbe_ = std::make_shared<amr::LoadBalancerEstimatorHybrid<PHARETypes>>(lb_info.mode,
                                                                               lbm_->getId());

    auto loadBalancer_db = std::make_shared<SAMRAI::tbox::MemoryDatabase>("LoadBalancerDB");
    loadBalancer_db->putDouble("flexible_load_tolerance", lb_info.tolerance);
    auto loadBalancer = std::make_shared<SAMRAI::mesh::CascadePartitioner>(
        SAMRAI::tbox::Dimension{dimension}, "LoadBalancer", loadBalancer_db);

    if (dict["simulation"]["AMR"]["refinement"].contains("tagging"))
    { // Load balancers break with refinement boxes - only tagging supported
        /*
          P=0000000:Program abort called in file ``/.../SAMRAI/xfer/RefineSchedule.cpp'' at line 369
          P=0000000:ERROR MESSAGE:
          P=0000000:RefineSchedule:RefineSchedule error: We are not currently
          P=0000000:supporting RefineSchedules with the source level finer
          P=0000000:than the destination level
        */
        lbm_->addLoadBalancerEstimator(maxMHDLevel_, maxLevelNumber_ - 1, std::move(lbe_));
        lbm_->setLoadBalancer(loadBalancer);
    }

    auto lbm_id = lbm_->getId(); // moved on next line
    multiphysInteg_->setLoadBalancerManager(std::move(lbm_));

    startTime_ = restart_time(dict);

    integrator_
        = std::make_unique<Integrator>(dict, hierarchy_, multiphysInteg_, multiphysInteg_,
                                       loadBalancer, startTime_, finalTime_, lb_info, lbm_id);

    timeStamper = core::TimeStamperFactory::create(dict["simulation"]);

    if (dict["simulation"].contains("diagnostics"))
        diagnostics_init(dict["simulation"]["diagnostics"], *hybridModel_);
}


template<auto opts>
void Simulator<opts>::mhd_init(initializer::PHAREDict const& dict)
{
    mhdModel_ = std::make_shared<MHDModel>(dict["simulation"], mhd_resman_ptr);
    mhd_resman_ptr->registerResources(mhdModel_->state);

    // we register the mhd model for all possible levels in the hierarchy
    // since for now it is the only model available, same for the solver
    multiphysInteg_->registerModel(0, maxMHDLevel_ - 1, mhdModel_);

    multiphysInteg_->registerAndInitSolver(0, maxMHDLevel_ - 1,
                                           std::make_unique<SolverMHD>(dict["simulation"]["algo"]));

    multiphysInteg_->registerAndSetupMessengers(messengerFactory_);

    if (dict["simulation"]["AMR"]["refinement"].contains("tagging"))
    {
        if (dict["simulation"]["AMR"]["refinement"]["tagging"]["method"].template to<std::string>()
            != "none")
        {
            auto mhdTagger_ = amr::TaggerFactory<MHDModel>::make(
                dict["simulation"]["AMR"]["refinement"]["tagging"]);
            multiphysInteg_->registerTagger(0, maxMHDLevel_ - 1, std::move(mhdTagger_));
        }
    }

    amr::LoadBalancerDetails lb_info
        = amr::LoadBalancerDetails::FROM(dict["simulation"]["AMR"]["loadbalancing"]);

    auto lbm_ = std::make_unique<amr::LoadBalancerManager<opts.dimension>>(dict);
    auto lbe_ = std::make_shared<amr::LoadBalancerEstimatorMHD<PHARETypes>>(lbm_->getId());

    auto loadBalancer_db = std::make_shared<SAMRAI::tbox::MemoryDatabase>("LoadBalancerDB");
    loadBalancer_db->putDouble("flexible_load_tolerance", lb_info.tolerance);
    auto loadBalancer = std::make_shared<SAMRAI::mesh::CascadePartitioner>(
        SAMRAI::tbox::Dimension{dimension}, "LoadBalancer", loadBalancer_db);

    if (dict["simulation"]["AMR"]["refinement"].contains("tagging"))
    { // Load balancers break with refinement boxes - only tagging supported
        /*
          P=0000000:Program abort called in file ``/.../SAMRAI/xfer/RefineSchedule.cpp'' at line 369
          P=0000000:ERROR MESSAGE:
          P=0000000:RefineSchedule:RefineSchedule error: We are not currently
          P=0000000:supporting RefineSchedules with the source level finer
          P=0000000:than the destination level
        */
        lbm_->addLoadBalancerEstimator(0, maxMHDLevel_ - 1, std::move(lbe_));
        lbm_->setLoadBalancer(loadBalancer);
    }

    auto lbm_id = lbm_->getId(); // moved on next line
    multiphysInteg_->setLoadBalancerManager(std::move(lbm_));

    startTime_ = restart_time(dict);

    integrator_
        = std::make_unique<Integrator>(dict, hierarchy_, multiphysInteg_, multiphysInteg_,
                                       loadBalancer, startTime_, finalTime_, lb_info, lbm_id);

    timeStamper = core::TimeStamperFactory::create(dict["simulation"]);

    if (dict["simulation"].contains("diagnostics"))
        diagnostics_init(dict["simulation"]["diagnostics"], *mhdModel_);
}



template<auto opts>
Simulator<opts>::Simulator(PHARE::initializer::PHAREDict const& dict,
                           std::shared_ptr<PHARE::amr::Hierarchy> const& hierarchy)
    : coutbuf{logging(log_out)}
    , hierarchy_{hierarchy}
    , modelNames_{dict["simulation"]["models"].template to<std::vector<std::string>>()}
    , descriptors_{PHARE::amr::makeDescriptors(modelNames_)}
    , messengerFactory_{descriptors_}
    , maxLevelNumber_{dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>()}
    , maxMHDLevel_{dict["simulation"]["AMR"]["max_mhd_level"].template to<int>()}
    , dt_{dict["simulation"]["time_step"].template to<double>()}
    , timeStepNbr_{dict["simulation"]["time_step_nbr"].template to<int>()}
    , finalTime_{dt_ * timeStepNbr_}
    , functors_{functors_setup(dict)}
    , multiphysInteg_{std::make_shared<MultiPhysicsIntegrator>(dict["simulation"], functors_)}
{
    currentTime_ = restart_time(dict);
    finalTime_ += currentTime_;


    // we would need a different restart manager for mhd and hybrid if both models are used

    if (find_model("HybridModel"))
    {
        hyb_resman_ptr = std::make_shared<HybridResourceManager_t>();
        hybrid_init(dict);
        if (dict["simulation"].contains("restarts"))
            rMan = restarts::RestartsManagerResolver::make_unique(*hierarchy_, *hyb_resman_ptr,
                                                                  dict["simulation"]["restarts"]);
    }

    if (find_model("MHDModel"))
    {
        mhd_resman_ptr = std::make_shared<MHDResourceManager_t>();
        mhd_init(dict);
        if (dict["simulation"].contains("restarts"))
            rMan = restarts::RestartsManagerResolver::make_unique(*hierarchy_, *mhd_resman_ptr,
                                                                  dict["simulation"]["restarts"]);
    }

    if (!hyb_resman_ptr and !mhd_resman_ptr)
        throw std::runtime_error("unsupported model");

    amr::ResourcesManagerGlobals::registerForRestarts();
}



template<auto opts>
std::string Simulator<opts>::to_str()
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




template<auto opts>
void Simulator<opts>::initialize()
{
    PHARE_LOG_SCOPE(1, "Simulator::initialize");

    std::optional<std::string> error = std::nullopt;

    try
    {
        if (isInitialized)
            throw std::runtime_error("cannot initialize  - simulator already isInitialized");

        if (integrator_ == nullptr)
            throw std::runtime_error("Error - Simulator has no integrator");

        integrator_->initialize();
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

    if (core::mpi::any(core::Errors::instance().any()))
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
double Simulator<opts>::advance(double dt)
{
    PHARE_LOG_SCOPE(1, "Simulator::advance");

    double dt_new                    = 0;
    std::optional<std::string> error = std::nullopt;

    try
    {
        if (!integrator_)
            throw std::runtime_error("Error - no valid integrator in the simulator");

        dt_new       = integrator_->advance(dt);
        currentTime_ = startTime_ + ((*timeStamper) += dt);
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

    if (core::mpi::any(core::Errors::instance().any()))
    {
        this->dMan.reset(); // closes/flushes hdf5 files
        if (error)
            throw std::runtime_error(*error);
        throw std::runtime_error("forcing error");
    }

    return dt_new;
}




template<auto opts>
auto Simulator<opts>::find_model(std::string name)
{
    if (modelNames_.empty())
        throw std::runtime_error("Simulator: No models found!");
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
            SimOpts<> constexpr static opts{
                .dimension = d, .interp_order = io, .nbRefinedPart = nb};
            return std::make_unique<Simulator<opts>>(theDict, hierarchy_);
        }
        else
        {
            return nullptr;
        }
    }
};


template<typename Simulator>
std::unique_ptr<Simulator> makeSimulator(std::shared_ptr<amr::Hierarchy> const& hierarchy)
{
    return std::make_unique<Simulator>(initializer::PHAREDictHandler::INSTANCE().dict(), hierarchy);
}



} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_H*/
