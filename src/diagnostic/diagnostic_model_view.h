#ifndef DIAGNOSTIC_MODEL_VIEW_H
#define DIAGNOSTIC_MODEL_VIEW_H

#include "core/utilities/mpi_utils.h"
#include "solver/physical_models/hybrid_model.h"
#include "cppdict/include/dict.hpp"

namespace PHARE::diagnostic
{
// Generic Template declaration, to override per Concrete model type
class IModelView
{
public:
    inline virtual ~IModelView();
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
    using GridLayout = typename Model::gridlayout_type;
    using PatchProperties
        = cppdict::Dict<float, double, std::size_t, std::vector<int>, std::vector<std::uint32_t>,
                        std::vector<double>, std::string>;

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

    auto domainBox() const { return hierarchy_.domainBox(); }

    auto origin() const { return hierarchy_.origin(); }

    auto cellWidth() const { return hierarchy_.cellWidth(); }

    std::string getLayoutTypeString() const { return std::string{GridLayout::implT::type}; }

    static auto getPatchProperties(GridLayout const& grid)
    {
        PatchProperties dict;
        dict["origin"]   = grid.origin().toVector();
        dict["nbrCells"] = core::Point<std::uint32_t, dimension>{grid.nbrCells()}.toVector();
        dict["lower"]    = grid.AMRBox().lower.toVector();
        dict["upper"]    = grid.AMRBox().upper.toVector();
        dict["mpi_rank"] = static_cast<std::size_t>(core::mpi::rank());
        return dict;
    }


    static auto getEmptyPatchProperties()
    {
        PatchProperties dict;
        dict["origin"]   = std::vector<double>{};
        dict["nbrCells"] = std::vector<std::uint32_t>{};
        dict["lower"]    = std::vector<int>{};
        dict["upper"]    = std::vector<int>{};
        dict["mpi_rank"] = std::size_t{0};
        return dict;
    }


protected:
    Model& model_;
    Hierarchy& hierarchy_;
};



} // namespace PHARE::diagnostic



#endif // DIAGNOSTIC_MODEL_VIEW_H
