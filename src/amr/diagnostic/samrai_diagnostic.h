#ifndef PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H
#define PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H

#include "types/amr_types.h"
#include "diagnostic_manager.h"
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/Patch.h>

namespace PHARE
{
// generic subclass of model specialized superclass
template<typename Model>
class SamraiDiagnosticModelView : public DiagnosticModelView<Model, typename Model::type_list>
{
public:
    using Super      = DiagnosticModelView<Model, typename Model::type_list>;
    using ResMan     = typename Model::resources_manager_type;
    using GridLayout = typename Model::gridLayout_type;
    using Guard      = amr::ResourcesGuard<ResMan, Model>;
    using Hierarchy  = amr::SAMRAI_Types::hierarchy_t;
    using Patch      = amr::SAMRAI_Types::patch_t;
    using Super::model_;
    static constexpr auto dimension = Model::dimension;

    SamraiDiagnosticModelView(Hierarchy& hierarchy, Model& model)
        : Super{model}
        , hierarchy_{hierarchy}
    {
    }

    auto guardedGrid(Patch& patch) { return GuardedGrid{patch, Super::model_}; }


    template<typename Action, typename... Args>
    void visitLevel(int level, Action&& action)
    {
        amr::visitLevel<GridLayout>(*hierarchy_.getPatchLevel(level), *model_.resourcesManager,
                                    action, model_);
    }

protected:
    struct GuardedGrid
    {
        using Guard = typename SamraiDiagnosticModelView<Model>::Guard;

        GuardedGrid(Patch& patch, Model& model)
            : guard_{model.resourcesManager->setOnPatch(patch, model)}
            , grid_{PHARE::amr::layoutFromPatch<GridLayout>(patch)}
        {
        }

        operator GridLayout&() { return grid_; }

        Guard guard_;
        GridLayout grid_;
    };


private:
    Hierarchy& hierarchy_;

    SamraiDiagnosticModelView(const SamraiDiagnosticModelView&)             = delete;
    SamraiDiagnosticModelView(const SamraiDiagnosticModelView&&)            = delete;
    SamraiDiagnosticModelView& operator&(const SamraiDiagnosticModelView&)  = delete;
    SamraiDiagnosticModelView& operator&(const SamraiDiagnosticModelView&&) = delete;
};


} /*namespace PHARE*/

#endif /*PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H*/
