#ifndef DIAGNOSTIC_MODEL_VIEW_HPP
#define DIAGNOSTIC_MODEL_VIEW_HPP

#include "core/def.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "cppdict/include/dict.hpp"
#include <type_traits>

namespace PHARE::diagnostic
{
// Generic Template declaration, to override per Concrete model type
class IModelView
{
public:
    inline virtual ~IModelView();
};
IModelView::~IModelView() {}

template<typename Hierarchy, typename Model>
class BaseModelView : public IModelView
{
public:
    using GridLayout = typename Model::gridlayout_type;
    using VecField   = typename Model::vecfield_type;
    using PatchProperties
        = cppdict::Dict<float, double, std::size_t, std::vector<int>, std::vector<std::uint32_t>,
                        std::vector<double>, std::vector<std::size_t>, std::string,
                        std::vector<std::string>>;

    BaseModelView(Hierarchy& hierarchy, Model& model)
        : model_{model}
        , hierarchy_{hierarchy}
    {
    }

    template<typename Action>
    void visitHierarchy(Action&& action, int minLevel = 0, int maxLevel = 0)
    {
        PHARE::amr::visitHierarchy<GridLayout>(hierarchy_, *model_.resourcesManager,
                                               std::forward<Action>(action), minLevel, maxLevel,
                                               model_);
    }

    NO_DISCARD auto boundaryConditions() const { return hierarchy_.boundaryConditions(); }
    NO_DISCARD auto domainBox() const { return hierarchy_.domainBox(); }
    NO_DISCARD auto origin() const { return hierarchy_.origin(); }
    NO_DISCARD auto cellWidth() const { return hierarchy_.cellWidth(); }

    NO_DISCARD std::string getLayoutTypeString() const
    {
        return std::string{GridLayout::implT::type};
    }

    NO_DISCARD auto getPatchProperties(std::string patchID, GridLayout const& grid) const
    {
        PatchProperties dict;
        dict["origin"]   = grid.origin().toVector();
        dict["nbrCells"] = core::Point<std::uint32_t, Model::dimension>{grid.nbrCells()}.toVector();
        dict["lower"]    = grid.AMRBox().lower.toVector();
        dict["upper"]    = grid.AMRBox().upper.toVector();
        dict["mpi_rank"] = static_cast<std::size_t>(core::mpi::rank());
        return dict;
    }

    NO_DISCARD static auto getEmptyPatchProperties(PatchProperties dict = {})
    {
        dict["origin"]   = std::vector<double>{};
        dict["nbrCells"] = std::vector<std::uint32_t>{};
        dict["lower"]    = std::vector<int>{};
        dict["upper"]    = std::vector<int>{};
        dict["mpi_rank"] = std::size_t{0};
        return dict;
    }

    NO_DISCARD bool hasTagsVectorFor(int ilevel, std::string patch_id) const
    {
        auto key = std::to_string(ilevel) + "_" + patch_id;
        return model_.tags.count(key);
    }

    NO_DISCARD auto& getTagsVectorFor(int ilevel, std::string patch_id) const
    {
        auto key = std::to_string(ilevel) + "_" + patch_id;
        return model_.tags.at(key);
    }

protected:
    Model& model_;
    Hierarchy& hierarchy_;
};

// we need to check the model type or the template specialization from the typelist will fail to
// compile for the other model.
template<typename Model>
struct is_hybrid_model : std::false_type
{
};

template<typename... Args>
struct is_hybrid_model<solver::HybridModel<Args...>> : std::true_type
{
};

template<typename Model>
struct is_mhd_model : std::false_type
{
};

template<typename... Args>
struct is_mhd_model<solver::MHDModel<Args...>> : std::true_type
{
};

template<typename HybridModel>
bool constexpr is_valid_hybrid_model = []() {
    if constexpr (is_hybrid_model<HybridModel>::value)
    {
        return std::is_same_v<solver::type_list_to_hybrid_model_t<typename HybridModel::type_list>,
                              HybridModel>;
    }
    else
    {
        return false;
    }
}();

template<typename MHDModel>
bool constexpr is_valid_mhd_model = []() {
    if constexpr (is_mhd_model<MHDModel>::value)
    {
        return std::is_same_v<solver::type_list_to_mhd_model_t<typename MHDModel::type_list>,
                              MHDModel>;
    }
    else
    {
        return false;
    }
}();

struct HybridIdentifier
{
};
struct MHDIdentifier
{
};

template<typename Model>
using ModelIdentifier
    = std::conditional_t<is_valid_hybrid_model<Model>, HybridIdentifier,
                         std::conditional_t<is_valid_mhd_model<Model>, MHDIdentifier, void>>;

template<typename Hierarchy, typename Model, typename = ModelIdentifier<Model>>
class ModelView
{
    static_assert(is_valid_hybrid_model<Model> || is_valid_mhd_model<Model>,
                  "Model must be either a valid HybridModel or MHDModel.");
};

template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, HybridIdentifier> : public BaseModelView<Hierarchy, Model>
{
    using VecField = typename Model::vecfield_type;

public:
    using Identifier = HybridIdentifier;

    using BaseModelView<Hierarchy, Model>::BaseModelView;

    NO_DISCARD VecField& getB() const { return this->model_.state.electromag.B; }

    NO_DISCARD VecField& getE() const { return this->model_.state.electromag.E; }

    NO_DISCARD auto& getIons() const { return this->model_.state.ions; }
};

template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, MHDIdentifier> : public BaseModelView<Hierarchy, Model>
{
    using Field    = typename Model::field_type;
    using VecField = typename Model::vecfield_type;

public:
    using Identifier = MHDIdentifier;

    using BaseModelView<Hierarchy, Model>::BaseModelView;

    NO_DISCARD Field& getRho() const { return this->model_.state.rho; }

    NO_DISCARD VecField& getV() const { return this->model_.state.V; }

    NO_DISCARD VecField& getB() const { return this->model_.state.B; }

    NO_DISCARD Field& getP() const { return this->model_.state.P; }

    NO_DISCARD VecField& getRhoV() const { return this->model_.state.rhoV; }

    NO_DISCARD Field& getEtot() const { return this->model_.state.Etot; }

    // What about E, J, etc. ? Those 2 in particular are physically relevant temporaries.
};


} // namespace PHARE::diagnostic



#endif // DIAGNOSTIC_MODEL_VIEW_HPP
