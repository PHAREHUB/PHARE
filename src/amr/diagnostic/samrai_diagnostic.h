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
<<<<<<< HEAD
    using Super = DiagnosticModelView<Model, typename Model::type_list>;
    using Grid  = typename Model::gridLayout_type;
    using Guard = typename Super::Guard;
    using Super::getResources;
=======
    using Super      = DiagnosticModelView<Model, typename Model::type_list>;
    using ResMan     = typename Model::resources_manager_type;
    using GridLayout = typename Model::gridLayout_type;
    using Guard      = amr::ResourcesGuard<ResMan, Model>;
    using Patch      = ::SAMRAI::hier::Patch;
>>>>>>> wip diag pr
    using Super::model_;
    using typename Super::ResMan;
    using typename Super::Resources;
    using Patch = ::SAMRAI::hier::Patch;

<<<<<<< HEAD
    SamraiModelDiagnosticView(Model& model)
        : Super(model)
=======
    SamraiDiagnosticModelView(Model& model)
        : Super{model}
>>>>>>> wip diag pr
    {
    }

    auto guardedGrid(Patch& patch)
    {
        return std::apply(
            [&patch, this](auto&... args) { return GuardedGrid(patch, model_, args...); },
            getResources());
    }

protected:
    struct GuardedGrid
    {
        using Guard = typename SamraiDiagnosticModelView<Model>::Guard;

<<<<<<< HEAD
        template<typename... Args>
        GuardedGrid(Patch& patch, Model& model, Args&... args)
            : guard_(model.resourcesManager->setOnPatch(patch, args...))
            , grid_(PHARE::amr::layoutFromPatch<Grid>(patch))
=======
        GuardedGrid(Patch& patch, Model& model)
            : guard_{model.resourcesManager->setOnPatch(patch, model)}
            , grid_{PHARE::amr::layoutFromPatch<GridLayout>(patch)}
>>>>>>> wip diag pr
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
