#ifndef PHARE_SOLVER_PPC_STEPS_HPP
#define PHARE_SOLVER_PPC_STEPS_HPP

// SolverPPC predictor1_/predictor2_/corrector_ (Faraday, Ampere, Ohm), kept out of
// solver_ppc.hpp so one translation unit does not instantiate the whole solver at once
// include this wherever SolverPPC is instantiated, or emit the steps once via SolverPPCSteps

#include "amr/solvers/solver_ppc.hpp"

#include <array>

namespace PHARE::solver
{

template<typename Solver>
struct SolverPPCSteps
{
    // forces the steps to be emitted in the calling translation unit
    static void emit()
    {
        [[gnu::used]] static constexpr std::array steps{&Solver::predictor1_, &Solver::predictor2_,
                                                        &Solver::corrector_};
    }
};


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor1_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_");

    TimeSetter setTime{level, model, newTime};

    Faraday_t faraday{level, model};
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.faraday");
        auto const dt = newTime - currentTime;
        faraday(model.state.electromag.B, model.state.electromag.E, electromagPred_.B, dt);
        setTime(electromagPred_.B);
    }
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.faraday::schedules");
        fromCoarser.fillMagneticGhosts(electromagPred_.B, level, newTime);
    }

    Ampere_t ampere{level, model};
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ampere");
        ampere(electromagPred_.B, model.state.J);
        setTime(model.state.J);
    }


    Ohm_t ohm{ohm_info, level, model};
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ohm");
        update_electrons(level, model);
        ohm(electromagPred_.B, electromagPred_.E, model.state.electrons);
        setTime(electromagPred_.E);
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor2_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_");

    TimeSetter setTime{level, model, newTime};

    Faraday_t faraday{level, model};
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.faraday");
        auto dt = newTime - currentTime;
        faraday(model.state.electromag.B, electromagAvg_.E, electromagPred_.B, dt);
        setTime(electromagPred_.B);
    }
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.faraday::schedules");
        fromCoarser.fillMagneticGhosts(electromagPred_.B, level, newTime);
    }

    Ampere_t ampere{level, model};
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ampere");
        ampere(electromagPred_.B, model.state.J);
        setTime(model.state.J);
    }

    Ohm_t ohm{ohm_info, level, model};
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ohm");
        update_electrons(level, model);
        ohm(electromagPred_.B, electromagPred_.E, model.state.electrons);
        setTime(electromagPred_.E);
    }
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::corrector_(level_t& level, HybridModel& model,
                                                   Messenger& fromCoarser, double const currentTime,
                                                   double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::corrector_");

    auto levelNumber = level.getLevelNumber();
    TimeSetter setTime{level, model, newTime};

    auto& electromag = model.state.electromag;
    Faraday_t faraday{level, model};
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.faraday");
        auto dt = newTime - currentTime;
        faraday(electromag.B, electromagAvg_.E, electromag.B, dt);
        setTime(model.state.electromag.B);
    }
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.faraday::schedules");
        fromCoarser.fillMagneticGhosts(model.state.electromag.B, level, newTime);
    }

    Ampere_t ampere{level, model};
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ampere");
        ampere(electromag.B, model.state.J);
        setTime(model.state.J);
    }

    Ohm_t ohm{ohm_info, level, model};
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ohm");

        update_electrons(level, model);
        ohm(electromag.B, electromag.E, model.state.electrons);
        setTime(model.state.electromag.E);
    }
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ohm::schedules");
        fromCoarser.fillElectricGhosts(model.state.electromag.E, level, newTime);
    }
}




} // namespace PHARE::solver

#endif /*PHARE_SOLVER_PPC_STEPS_HPP*/
