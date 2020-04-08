#ifndef DIAGNOSTIC_MODEL_VIEW_H
#define DIAGNOSTIC_MODEL_VIEW_H

#include "solver/physical_models/hybrid_model.h"

namespace PHARE::diagnostic
{
// Generic Template declaration, to override per Concrete model type
class IModelView
{
public:
    virtual ~IModelView();
};
IModelView::~IModelView() {}



template<typename Model>
bool constexpr is_hybrid_model
    = std::is_same_v<solver::type_list_to_hybrid_model_t<typename Model::type_list>, Model>;




template<typename Hierarchy, typename Model, std::enable_if_t<is_hybrid_model<Model>, int> = 0>
class ModelView : public IModelView
{
    using VecField                  = typename Model::vecfield_type;
    using ResMan                    = typename Model::resources_manager_type;
    static constexpr auto dimension = Model::dimension;


public:
    using GridLayout      = typename Model::gridLayout_type;
    using PatchProperties = cppdict::Dict<float, double, size_t, std::string>;

    ModelView(Hierarchy& hierarchy, Model& model)
        : model_{model}
        , hierarchy_{hierarchy}
    {
    }

    std::vector<VecField*> getElectromagFields() const
    {
        return {&model_.state.electromag.B, &model_.state.electromag.E};
    }

    auto& getIons() const { return model_.state.ions; }


    template<typename Action, typename... Args>
    void visitHierarchy(Action&& action, int minLevel = 0, int maxLevel = 0)
    {
        auto& resMan = *model_.resourcesManager;
        PHARE::amr::visitHierarchy<GridLayout>(hierarchy_, resMan, std::forward<Action>(action),
                                               minLevel, maxLevel, model_);
    }

    std::string getLayoutTypeString() const { return std::string{GridLayout::implT::type}; }

    std::string getMeshSize() const { return hierarchy_.meshSize(); }

    std::string getCells() const { return hierarchy_.cells(); }

    std::string getOrigin() const { return hierarchy_.origin(); }

    static auto getPatchProperties(GridLayout const& grid)
    {
        PatchProperties dict;
        dict["origin"]   = grid.origin().str();
        dict["nbrCells"] = core::Point<uint32_t, dimension>{grid.nbrCells()}.str();
        return dict;
    }


    static auto getEmptyPatchProperties()
    {
        PatchProperties dict;
        dict["origin"]   = std::string{""};
        dict["nbrCells"] = std::string{""};
        return dict;
    }


protected:
    Model& model_;
    Hierarchy& hierarchy_;
};




} // namespace PHARE::diagnostic




#endif // DIAGNOSTIC_MODEL_VIEW_H
