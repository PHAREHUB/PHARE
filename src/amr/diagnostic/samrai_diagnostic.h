#ifndef PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H
#define PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H

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
    using Patch      = ::SAMRAI::hier::Patch;
    using Super::model_;

    SamraiDiagnosticModelView(Model& model)
        : Super{model}
    {
    }

    auto guardedGrid(Patch& patch) { return GuardedGrid(patch, Super::model_); }

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
    SamraiDiagnosticModelView(const SamraiDiagnosticModelView&)             = delete;
    SamraiDiagnosticModelView(const SamraiDiagnosticModelView&&)            = delete;
    SamraiDiagnosticModelView& operator&(const SamraiDiagnosticModelView&)  = delete;
    SamraiDiagnosticModelView& operator&(const SamraiDiagnosticModelView&&) = delete;
};


template<typename ModelView>
class SamraiDiagnostic
{
public:
    using Hierarchy  = SAMRAI::hier::PatchHierarchy;
    using PatchLevel = std::shared_ptr<SAMRAI::hier::PatchLevel>;

    auto& modelView() { return modelView_; }

protected:
    SamraiDiagnostic(Hierarchy& hierarchy, ModelView& model)
        : modelView_{model}
        , hierarchy_{hierarchy}
    {
    }

    SamraiDiagnosticModelView<ModelView> modelView_;
    Hierarchy& hierarchy_;

private:
    SamraiDiagnostic(const SamraiDiagnostic&)             = delete;
    SamraiDiagnostic(const SamraiDiagnostic&&)            = delete;
    SamraiDiagnostic& operator&(const SamraiDiagnostic&)  = delete;
    SamraiDiagnostic& operator&(const SamraiDiagnostic&&) = delete;
};

} /*namespace PHARE*/

#endif /*PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H*/
