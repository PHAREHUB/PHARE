

#ifndef PHARE_SOLVER_PPC_H
#define PHARE_SOLVER_PPC_H

#include <SAMRAI/hier/Patch.h>


#include "evolution/solvers/solver.h"
#include "evolution/transactions/hybrid_transaction.h"
#include "evolution/transactions/hybrid_transaction_info.h"

namespace PHARE
{
template<typename HybridModel>
class SolverPPC : public ISolver
{
private:
    using Electromag = decltype(std::declval<HybridModel>().state.electromag);
    using IonsT      = decltype(std::declval<HybridModel>().state.ions);
    using VecFieldT  = decltype(std::declval<HybridModel>().state.electromag.E);

    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};


public:
    explicit SolverPPC()
        : ISolver{"PPC"}
    {
    }

    virtual ~SolverPPC() = default;

    virtual std::string modelName() const override { return HybridModel::model_name; }


    virtual void fillTransactionInfo(std::unique_ptr<ITransactionInfo> const& info) const override
    {
        auto& modelInfo = dynamic_cast<HybridTransactionInfo&>(*info);

        auto const& Epred = electromagPred_.E;
        auto const& Bpred = electromagPred_.B;

        auto magneticComponentNames = extractNames(Bpred);
        auto electricComponentNames = extractNames(Epred);

        modelInfo.ghostElectric.push_back({Epred.name(), electricComponentNames[0],
                                           electricComponentNames[1], electricComponentNames[2]});
        modelInfo.ghostMagnetic.push_back({Bpred.name(), magneticComponentNames[0],
                                           magneticComponentNames[1], magneticComponentNames[2]});

        // here fill modelInfo.initElectric and initMagnetic if some solver
    }


    virtual void registerResources(PhysicalModel& model) override
    {
        auto& hmodel = dynamic_cast<HybridModel&>(model);
        hmodel.resourcesManager->registerResources(electromagPred_.E);
        hmodel.resourcesManager->registerResources(electromagPred_.B);
        hmodel.resourcesManager->registerResources(electromagAvg_.E);
        hmodel.resourcesManager->registerResources(electromagAvg_.B);
    }

    virtual void allocate(PhysicalModel& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const override
    {
        auto& hmodel = dynamic_cast<HybridModel&>(model);
        hmodel.resourcesManager->allocate(electromagPred_.E, patch, allocateTime);
        hmodel.resourcesManager->allocate(electromagPred_.B, patch, allocateTime);
        hmodel.resourcesManager->allocate(electromagAvg_.E, patch, allocateTime);
        hmodel.resourcesManager->allocate(electromagAvg_.B, patch, allocateTime);
    }
    /*
        template<typename HybridTransaction>
        void initTransaction(HybridTransaction& hybridTransaction)
        {
            bool constexpr withTemporal{true};
            hybridTransaction.registerMagneticIn(electromagPred_.B,
       BooleanSelector<withTemporal>{}); hybridTransaction.registerElectricIn(electromagPred_.E,
       BooleanSelector<withTemporal>{});

            hybridTransaction.registerMagneticIn(electromagAvg_.B, BooleanSelector<withTemporal>{});
            hybridTransaction.registerElectricIn(electromagAvg_.E, BooleanSelector<withTemporal>{});

            hybridTransaction.registerMagneticSync(model_.electromag.B);
            hybridTransaction.registerElectricSync(model_.electromag.E);
        }*/


    virtual void advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                              int const levelNumber, PhysicalModel& model,
                              ITransaction& fromCoarserTransaction, const double currentTime,
                              const double newTime) override
    {
        // bool constexpr withTemporal{true};

        auto& hybridModel = dynamic_cast<HybridModel&>(model);
        auto& hybridState = hybridModel.state;
        auto& fromCoarser = dynamic_cast<HybridTransaction<HybridModel>&>(fromCoarserTransaction);


        /*
         * PREDICTOR 1
         */

        // loop on patches
        // |
        // -> faraday E , B, Bpred
        VecFieldT& Bpred = electromagPred_.B;
        fromCoarser.fillMagneticGhosts(Bpred, levelNumber, newTime);


        // loop on patches
        // |
        // -> ampere Bpred, Jtot on interior + ghost
        // -> ohm Bpred, ions.rho, Vepred1, PePred1, Jtot, Epred

        VecFieldT& Epred = electromagPred_.E;
        fromCoarser.fillElectricGhosts(Epred, levelNumber, newTime);
        // fromCoarser.getElectric(electromagPred_.E, fillTime, BooleanSelector<withTemporal>{},
        //                        FillTypeSelector<FillType::GhostRegion>{});

        // loop on patches
        // |
        // -> timeAverage E, Epred, Eavg
        // -> timeAverage B, Bpred, Bavg

        // fill PRA and ghostsin purple region so that some of these
        // particles may eventually enter the domain after being pushed
        fromCoarser.fillIonGhostParticles(hybridState.ions, levelNumber, newTime);

        // move all ions in:
        //  - the domain
        //  - the intra level ghost region
        //  - coarse to fine ghost region
        // accumulate those that are within the domain after being pushed (before pivot)
        // see particle design code

        // now some of the nodes close to the boundary of the patches (and level)
        // are incomplete because they may have recieved contributions from particles outside
        // their domain that would have entered the purple region during the push
        // therefore we need to re-fill the purple region and accumulate that density
        // this is done by calling a transaction to fill the moments.

        fromCoarser.fillIonMomentGhosts(hybridState.ions, levelNumber, newTime);




        // move ions -> needs the fromCoarser and toFiner but will also needs some real boundary
        // condition
        // needs predictor1 boolean

        /*
         * PREDICTOR 2
         */

        // loop on patches
        // |
        // -> faraday Eavg , B, Bpred

        fromCoarser.fillMagneticGhosts(Bpred, levelNumber, newTime);

        // loop on patches
        // |
        // -> ampere Bpred, Jtot on interior + ghost
        // -> ohm Bpred, ions.rho, Vepred2, PePred2, Jtot, Epred

        fromCoarser.fillElectricGhosts(Epred, levelNumber, newTime);


        // loop on patches
        // |
        // -> timeAverage E, Epred, Eavg
        // -> timeAverage B, Bpred, Bavg


        // fill PRA and ghostsin purple region so that some of these
        // particles may eventually enter the domain after being pushed
        fromCoarser.fillIonGhostParticles(hybridState.ions, levelNumber, newTime);


        // same as for predictor 1, except that we will updates the ions
        // populations

        fromCoarser.fillIonMomentGhosts(hybridState.ions, levelNumber, newTime);


        // needs predictor2 boolean


        /*
         * CORRECTOR
         */

        // loop on patches
        // |
        // -> faraday Eavg , B, B

        VecFieldT& B = hybridState.electromag.B;
        fromCoarser.fillMagneticGhosts(B, levelNumber, newTime);

        // loop on patches
        // |
        // -> ampere B, Jtot on interior + ghost
        // -> ohm Bpred, ions.rho, Vecorr, Pecorr, Jtot, E


        VecFieldT& E = hybridState.electromag.E;
        fromCoarser.fillElectricGhosts(E, levelNumber, newTime);


        // double newTime = 0.0;
        // return newTime;
    }
    /*
    template<typename HybridTransaction>
    void syncLevel(HybridTransaction& toCoarser)
    {
        toCoarser.syncMagnetic(model_.electromag.B);
        toCoarser.syncElectric(model_.electromag.E);
    }*/
};

} // namespace PHARE




/*
template<bool>
struct BooleanSelector
{
};

template<typename Enum, Enum>
class EnumSelector
{
};

enum class FillType {
    SameLevel,
    GhostRegionOnSameLevel,
    LevelBorderOnly,
    GhostRegion,
    EraseDestination
};



template<FillType fillValue>
using FillTypeSelector = EnumSelector<FillType, fillValue>;
*/


#endif
