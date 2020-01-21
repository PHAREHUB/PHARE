

#ifndef PHARE_SOLVER_PPC_H
#define PHARE_SOLVER_PPC_H

#include <SAMRAI/hier/Patch.h>

#include "core/data/particles/particle_array.h"
#include "initializer/data_provider.h"
#include "amr/messengers/hybrid_messenger.h"
#include "amr/messengers/hybrid_messenger_info.h"
#include "amr/resources_manager/amr_utils.h"
#include "solver/solvers/solver.h"
#include "core/numerics/pusher/pusher.h"
#include "core/numerics/pusher/pusher_factory.h"
#include "core/numerics/interpolator/interpolator.h"
#include "core/utilities/particle_selector/particle_selector.h"
#include "core/numerics/boundary_condition/boundary_condition.h"


namespace PHARE
{
namespace solver
{
    template<typename HybridModel, typename AMR_Types>
    class SolverPPC : public ISolver<AMR_Types>
    {
    private:
        static constexpr auto dimension    = HybridModel::dimension;
        static constexpr auto interp_order = HybridModel::gridLayout_type::interp_order;

        using Electromag = decltype(std::declval<HybridModel>().state.electromag);
        using Ions       = decltype(std::declval<HybridModel>().state.ions);
        using VecFieldT  = decltype(std::declval<HybridModel>().state.electromag.E);
        using GridLayout = typename HybridModel::gridLayout_type;



        Electromag electromagPred_{"EMPred"};
        Electromag electromagAvg_{"EMAvg"};


    public:
        using patch_t = typename AMR_Types::patch_t;
        using level_t = typename AMR_Types::level_t;

        explicit SolverPPC(PHARE::initializer::PHAREDict dict)
            : ISolver<AMR_Types>{"PPC"}
        {
        }

        virtual ~SolverPPC() = default;


        virtual std::string modelName() const override { return HybridModel::model_name; }

        virtual void
        fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


        virtual void registerResources(IPhysicalModel<AMR_Types>& model) override;


        virtual void allocate(IPhysicalModel<AMR_Types>& model, SAMRAI::hier::Patch& patch,
                              double const allocateTime) const override;



        virtual void advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                  int const levelNumber, IPhysicalModel<AMR_Types>& model,
                                  amr::IMessenger<IPhysicalModel<AMR_Types>>& fromCoarserMessenger,
                                  const double currentTime, const double newTime) override;



    private:
        enum class PredictorStep { predictor1, predictor2 };


        /*
        template<typename HybridMessenger>
        void syncLevel(HybridMessenger& toCoarser)
        {
            toCoarser.syncMagnetic(model_.electromag.B);
            toCoarser.syncElectric(model_.electromag.E);
        }*/


    }; // end solverPPC




    template<typename HybridModel, typename AMR_Types>
    void SolverPPC<HybridModel, AMR_Types>::registerResources(IPhysicalModel<AMR_Types>& model)
    {
        auto& hmodel = dynamic_cast<HybridModel&>(model);
        hmodel.resourcesManager->registerResources(electromagPred_);
        hmodel.resourcesManager->registerResources(electromagAvg_);
    }




    template<typename HybridModel, typename AMR_Types>
    void SolverPPC<HybridModel, AMR_Types>::allocate(IPhysicalModel<AMR_Types>& model,
                                                     SAMRAI::hier::Patch& patch,
                                                     double const allocateTime) const
    {
        auto& hmodel = dynamic_cast<HybridModel&>(model);
        hmodel.resourcesManager->allocate(electromagPred_, patch, allocateTime);
        hmodel.resourcesManager->allocate(electromagAvg_, patch, allocateTime);
    }




    template<typename HybridModel, typename AMR_Types>
    void SolverPPC<HybridModel, AMR_Types>::fillMessengerInfo(
        std::unique_ptr<amr::IMessengerInfo> const& info) const
    {
        auto& modelInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

        auto const& Epred = electromagPred_.E;
        auto const& Bpred = electromagPred_.B;

        modelInfo.ghostElectric.emplace_back(Epred);
        modelInfo.ghostMagnetic.emplace_back(Bpred);
    }




    template<typename HybridModel, typename AMR_Types>
    void SolverPPC<HybridModel, AMR_Types>::advanceLevel(
        std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy, int const levelNumber,
        IPhysicalModel<AMR_Types>& model,
        amr::IMessenger<IPhysicalModel<AMR_Types>>& fromCoarserMessenger, const double currentTime,
        const double newTime)
    {
        // bool constexpr withTemporal{true};

        auto& hybridModel = dynamic_cast<HybridModel&>(model);
        auto& hybridState = hybridModel.state;
        auto& fromCoarser
            = dynamic_cast<amr::HybridMessenger<HybridModel, IPhysicalModel<AMR_Types>>&>(
                fromCoarserMessenger);

        auto level = hierarchy->getPatchLevel(levelNumber);
        /*
         * PREDICTOR 1
         */

        // loop on patches
        // |
        // -> faraday E , B, Bpred
        // auto faraday = []() {};
        // amr::visitLevel(*level, *hybridModel.resourcesManager, faraday, hybridModel);

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


        // loop on patches
        // |
        // -> moveIons_(Bavg, Eavg, ions, PredictorStep::predictor1);


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

        auto coarserLevel = hierarchy->getPatchLevel(levelNumber - 1);
        fromCoarser.fillIonMomentGhosts(hybridState.ions, *level, currentTime, newTime);




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

        fromCoarser.fillIonMomentGhosts(hybridState.ions, *level, currentTime, newTime);


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




} // namespace solver


} // namespace PHARE

#endif
