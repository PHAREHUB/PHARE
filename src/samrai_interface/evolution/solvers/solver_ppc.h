

#ifndef PHARE_SOLVER_PPC_H
#define PHARE_SOLVER_PPC_H

#include <SAMRAI/hier/Patch.h>


#include "evolution/messengers/hybrid_messenger.h"
#include "evolution/messengers/hybrid_messenger_info.h"
#include "evolution/solvers/solver.h"

namespace PHARE
{
namespace amr_interface
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




        virtual void fillMessengerInfo(std::unique_ptr<IMessengerInfo> const& info) const override
        {
            auto& modelInfo = dynamic_cast<HybridMessengerInfo&>(*info);

            auto const& Epred = electromagPred_.E;
            auto const& Bpred = electromagPred_.B;

            modelInfo.ghostElectric.emplace_back(Epred);
            modelInfo.ghostMagnetic.emplace_back(Bpred);
        }




        virtual void registerResources(IPhysicalModel& model) override
        {
            auto& hmodel = dynamic_cast<HybridModel&>(model);
            hmodel.resourcesManager->registerResources(electromagPred_);
            hmodel.resourcesManager->registerResources(electromagAvg_);
        }




        virtual void allocate(IPhysicalModel& model, SAMRAI::hier::Patch& patch,
                              double const allocateTime) const override
        {
            auto& hmodel = dynamic_cast<HybridModel&>(model);
            hmodel.resourcesManager->allocate(electromagPred_, patch, allocateTime);
            hmodel.resourcesManager->allocate(electromagAvg_, patch, allocateTime);
        }




        virtual void advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                  int const levelNumber, IPhysicalModel& model,
                                  IMessenger& fromCoarserMessenger, const double currentTime,
                                  const double newTime) override
        {
            // bool constexpr withTemporal{true};

            auto& hybridModel = dynamic_cast<HybridModel&>(model);
            auto& hybridState = hybridModel.state;
            auto& fromCoarser = dynamic_cast<HybridMessenger<HybridModel>&>(fromCoarserMessenger);

            auto level = hierarchy->getPatchLevel(levelNumber);
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
            fromCoarser.fillIonGhostParticles(hybridState.ions, *level, newTime);

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
            // this is done by calling a messenger to fill the moments.

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
            fromCoarser.fillIonGhostParticles(hybridState.ions, *level, newTime);


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
        template<typename HybridMessenger>
        void syncLevel(HybridMessenger& toCoarser)
        {
            toCoarser.syncMagnetic(model_.electromag.B);
            toCoarser.syncElectric(model_.electromag.E);
        }*/
    };

} // namespace amr_interface


} // namespace PHARE

#endif
