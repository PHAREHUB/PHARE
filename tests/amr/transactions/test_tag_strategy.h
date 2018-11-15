#ifndef PHARE_TEST_TAG_STRATEGY_H
#define PHARE_TEST_TAG_STRATEGY_H


#include "data/field/field_data.h"
#include "data/grid/gridlayoutdefs.h"
#include "data/vecfield/vecfield_component.h"
#include "evolution/solvers/solver_ppc.h"
#include "evolution/transactions/hybrid_transaction.h"
#include "physical_models/hybrid_model.h"
#include "tools/amr_utils.h"

#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>

#include <map>
#include <string>



template<typename Field, typename GridLayout, typename Func>
void fillField(Field& field, GridLayout& layout, Func f)
{
    auto iStart = layout.physicalStartIndex(field, PHARE::Direction::X);
    auto iEnd   = layout.physicalEndIndex(field, PHARE::Direction::X);

    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        auto origin = layout.origin();
        auto x      = layout.fieldNodeCoordinates(field, origin, ix);
        field(ix)   = f(x[0]);
    }
}




template<typename HybridModel>
class TagStrategy : public SAMRAI::mesh::StandardTagAndInitStrategy
{
private:
    std::shared_ptr<HybridModel> model_;
    std::shared_ptr<PHARE::SolverPPC<HybridModel>> solver_;
    std::shared_ptr<PHARE::HybridTransaction<HybridModel>> transaction_;

public:
    explicit TagStrategy(std::shared_ptr<HybridModel> model,
                         std::shared_ptr<PHARE::SolverPPC<HybridModel>> solver,
                         std::shared_ptr<PHARE::HybridTransaction<HybridModel>> transaction)
        : model_{std::move(model)}
        , solver_{std::move(solver)}
        , transaction_{std::move(transaction)}
    {
        auto infoFromFiner   = transaction_->emptyInfoFromFiner();
        auto infoFromCoarser = transaction_->emptyInfoFromCoarser();

        model_->fillTransactionInfo(infoFromFiner);
        model_->fillTransactionInfo(infoFromCoarser);
        solver_->fillTransactionInfo(infoFromFiner);

        transaction_->registerQuantities(std::move(infoFromFiner), std::move(infoFromCoarser));
    }


    void initializeLevelData(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                             int const levelNumber, double const initDataTime,
                             bool const canBeRefined, bool const initialTime,
                             std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel
                             = std::shared_ptr<SAMRAI::hier::PatchLevel>(),
                             bool const allocateData = true) override
    {
        //


        auto level = hierarchy->getPatchLevel(levelNumber);

        if (allocateData)
        {
            for (auto patch : *level)
            {
                model_->allocate(*patch, initDataTime);
                solver_->allocate(*model_, *patch, initDataTime);
                transaction_->allocate(*patch, initDataTime);
            }
        }


        transaction_->registerLevel(hierarchy, levelNumber);


        if (oldLevel)
        {
            // in case of a regrid we need to make a bunch of temporary regriding schedules
            // using the init algorithms and actually perform the .fillData() for all of them
            transaction_->regrid(hierarchy, levelNumber, oldLevel, initDataTime);
        }


        else // we're creating a brand new finest level in the hierarchy
        {
            if (levelNumber >= 0)
            {
                // initializer.init(model);

                for (auto& patch : *level)
                {
                    auto _ = model_->resourcesManager->makeResourcesGuard(*patch,
                                                                          model_->state.electromag);

                    auto layout
                        = PHARE::layoutFromPatch<typename HybridModel::gridLayout_type>(*patch);

                    auto& Ex = model_->state.electromag.E.getComponent(PHARE::Component::X);
                    auto& Ey = model_->state.electromag.E.getComponent(PHARE::Component::Y);
                    auto& Ez = model_->state.electromag.E.getComponent(PHARE::Component::Z);

                    auto& Bx = model_->state.electromag.B.getComponent(PHARE::Component::X);
                    auto& By = model_->state.electromag.B.getComponent(PHARE::Component::Y);
                    auto& Bz = model_->state.electromag.B.getComponent(PHARE::Component::Z);


                    auto fillLevel
                        = [levelNumber](double) { return static_cast<double>(levelNumber); };


                    fillField(Ex, layout, fillEx);
                    fillField(Ey, layout, fillEy);
                    fillField(Ez, layout, fillEz);

                    fillField(Bx, layout, fillBx);
                    fillField(By, layout, fillBy);
                    fillField(Bz, layout, fillBz);
                }
            }

            else
            {
                transaction_->initLevel(levelNumber, initDataTime);
            }
        }
    }

    void resetHierarchyConfiguration(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                     int const coarsestLevel, int const finestLevel) override
    {
    }


    static double fillEx(double x) { return x; }
    static double fillEy(double x) { return 2 * x; }
    static double fillEz(double x) { return 3 * x; }

    static double fillBx(double x) { return 4 * x; }
    static double fillBy(double x) { return 5 * x; }
    static double fillBz(double x) { return 6 * x; }


    static double fillInt(double x) { return 1.; }

private:
};

#endif
