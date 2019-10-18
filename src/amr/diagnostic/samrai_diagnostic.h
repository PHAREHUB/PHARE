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
class SamraiModelDiagnosticView : public DiagnosticModelView<Model, typename Model::type_list>
{
public:
    using Super = DiagnosticModelView<Model, typename Model::type_list>;
    using Grid  = typename Model::gridLayout_type;
    using Guard = typename Super::Guard;
    using Super::getResources;
    using Super::model_;
    using typename Super::ResMan;
    using typename Super::Resources;
    using Patch = ::SAMRAI::hier::Patch;

    SamraiModelDiagnosticView(Model& model)
        : Super{model}
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
        using Guard = typename SamraiModelDiagnosticView<Model>::Guard;

        template<typename... Args>
        GuardedGrid(Patch& patch, Model& model, Args&... args)
            : guard_(model.resourcesManager->setOnPatch(patch, args...))
            , grid_(PHARE::amr::layoutFromPatch<Grid>(patch))
        {
        }

        operator Grid&() { return grid_; }

        Guard guard_;
        Grid grid_;
    };

private:
    SamraiModelDiagnosticView(const SamraiModelDiagnosticView&)             = delete;
    SamraiModelDiagnosticView(const SamraiModelDiagnosticView&&)            = delete;
    SamraiModelDiagnosticView& operator&(const SamraiModelDiagnosticView&)  = delete;
    SamraiModelDiagnosticView& operator&(const SamraiModelDiagnosticView&&) = delete;
};


template<typename Model>
class SamraiDiagnostic
{
public:
    using Hierarchy  = SAMRAI::hier::PatchHierarchy;
    using PatchLevel = std::shared_ptr<SAMRAI::hier::PatchLevel>;

    auto& modelView() { return modelView_; }

protected:
    SamraiDiagnostic(Hierarchy& hierarchy, Model& model)
        : modelView_(model)
        , hierarchy_(hierarchy)
    {
    }

    SamraiModelDiagnosticView<Model> modelView_;
    Hierarchy& hierarchy_;

private:
    SamraiDiagnostic(const SamraiDiagnostic&)             = delete;
    SamraiDiagnostic(const SamraiDiagnostic&&)            = delete;
    SamraiDiagnostic& operator&(const SamraiDiagnostic&)  = delete;
    SamraiDiagnostic& operator&(const SamraiDiagnostic&&) = delete;
};

} /*namespace PHARE*/

#endif /*PHARE_AMR_DIAGNOSTIC_SAMRAI_DIAGNOSTIC_H*/
