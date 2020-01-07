#ifndef PHARE_TEST_TAG_STRATEGY_H
#define PHARE_TEST_TAG_STRATEGY_H


#include "amr/types/amr_types.h"
#include "amr/data/field/field_data.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield_component.h"
#include "amr/messengers/hybrid_messenger.h"
#include "solver/physical_models/hybrid_model.h"
#include "amr/resources_manager/amr_utils.h"
#include "solver/solvers/solver_ppc.h"

#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>

#include <map>
#include <string>


using namespace PHARE::core;
using namespace PHARE::amr;
using namespace PHARE::solver;

template<typename Field, typename GridLayout, typename Func>
void fillField(Field& field, GridLayout& layout, Func f)
{
    auto iStart = layout.physicalStartIndex(field, Direction::X);
    auto iEnd   = layout.physicalEndIndex(field, Direction::X);

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
    std::shared_ptr<SolverPPC<HybridModel, SAMRAI_Types>> solver_;
    std::shared_ptr<HybridMessenger<HybridModel, IPhysicalModel<SAMRAI_Types>>> messenger_;

public:
    explicit TagStrategy(
        std::shared_ptr<HybridModel> model,
        std::shared_ptr<SolverPPC<HybridModel, SAMRAI_Types>> solver,
        std::shared_ptr<HybridMessenger<HybridModel, IPhysicalModel<SAMRAI_Types>>> messenger)
        : model_{std::move(model)}
        , solver_{std::move(solver)}
        , messenger_{std::move(messenger)}
    {
        auto infoFromFiner   = messenger_->emptyInfoFromFiner();
        auto infoFromCoarser = messenger_->emptyInfoFromCoarser();

        model_->fillMessengerInfo(infoFromFiner);
        model_->fillMessengerInfo(infoFromCoarser);
        solver_->fillMessengerInfo(infoFromFiner);

        messenger_->registerQuantities(std::move(infoFromFiner), std::move(infoFromCoarser));
    }


    void initializeLevelData(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                             int const levelNumber, double const initDataTime,
                             bool const canBeRefined, bool const initialTime,
                             std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel
                             = std::shared_ptr<SAMRAI::hier::PatchLevel>(),
                             bool const allocateData = true) override
    {
        auto level = hierarchy->getPatchLevel(levelNumber);
        if (allocateData)
        {
            for (auto patch : *level)
            {
                model_->allocate(*patch, initDataTime);
                solver_->allocate(*model_, *patch, initDataTime);
                messenger_->allocate(*patch, initDataTime);
            }
        }

        messenger_->registerLevel(hierarchy, levelNumber);

        if (oldLevel)
        {
            // in case of a regrid we need to make a bunch of temporary regriding schedules
            // using the init algorithms and actually perform the .fillData() for all of them
            messenger_->regrid(hierarchy, levelNumber, oldLevel, *model_, initDataTime);
        }


        else // we're creating a brand new finest level in the hierarchy
        {
            if (levelNumber == 0)
            {
                model_->initialize(*level);
                messenger_->fillRootGhosts(*model_, *level, initDataTime);
            }

            else
            {
                messenger_->initLevel(*model_, *level, initDataTime);
            }
        }

        if (levelNumber == 0)
        {
            messenger_->fillIonGhostParticles(model_->state.ions, *level, initDataTime);
        }
    }

    void resetHierarchyConfiguration(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                     int const coarsestLevel, int const finestLevel) override
    {
    }


    // it is important these functions are linear functions of x
    // we compare the values obtained on the refined level to the value returned
    // by these functions. This only works if they are linear functions of x
    // since the refinement operator is a linear interpolation
    static double fillEx(double x) { return x; }
    static double fillEy(double x) { return 2 * x; }
    static double fillEz(double x) { return 3 * x; }

    static double fillBx(double x) { return 4 * x; }
    static double fillBy(double x) { return 5 * x; }
    static double fillBz(double x) { return 6 * x; }

    static double fillInt([[maybe_unused]] double x) { return 1.; }

private:
};

#endif
