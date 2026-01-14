#ifndef DIAGNOSTIC_MODEL_VIEW_HPP
#define DIAGNOSTIC_MODEL_VIEW_HPP

#include "core/def.hpp"
#include "core/utilities/mpi_utils.hpp"

#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/messengers/field_sum_transaction.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"

#include "cppdict/include/dict.hpp"
#include <tuple>
#include <type_traits>

#include <SAMRAI/xfer/RefineAlgorithm.h>

#include <type_traits>
#include <utility>

namespace PHARE::diagnostic
{
// Generic Template declaration, to override per Concrete model type
class IModelView
{
public:
    inline virtual ~IModelView();
};
IModelView::~IModelView() {}


template<typename Derived, typename Hierarchy, typename Model>
class BaseModelView : public IModelView
{
public:
    using GridLayout        = Model::gridlayout_type;
    using VecField          = Model::vecfield_type;
    using GridLayoutT       = Model::gridlayout_type;
    using ResMan            = Model::resources_manager_type;
    using TensorFieldData_t = ResMan::template UserTensorField_t</*rank=*/2>::patch_data_type;
    static constexpr auto dimension = Model::dimension;

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
    void onLevels(Action&& action, int minlvl = 0, int maxlvl = 0)
    {
        for (int ilvl = minlvl; ilvl < hierarchy_.getNumberOfLevels() && ilvl <= maxlvl; ++ilvl)
            if (auto lvl = hierarchy_.getPatchLevel(ilvl))
                action(*lvl);
    }

    template<typename Action>
    void visitHierarchy(Action&& action, int minLevel = 0, int maxLevel = 0)
    {
        PHARE::amr::visitHierarchy<GridLayout>(hierarchy_, *model_.resourcesManager,
                                               std::forward<Action>(action), minLevel, maxLevel,
                                               model_, *this);
    }

    NO_DISCARD auto boundaryConditions() const { return hierarchy_.boundaryConditions(); }
    NO_DISCARD auto domainBox() const { return hierarchy_.domainBox(); }
    NO_DISCARD auto origin() const { return std::vector<double>(dimension, 0); }
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

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return derived().getCompileTimeResourcesViewList();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return derived().getCompileTimeResourcesViewList();
    }

protected:
    Model& model_;
    Hierarchy& hierarchy_;

private:
    Derived& derived() { return static_cast<Derived&>(*this); }
    Derived const& derived() const { return static_cast<Derived const&>(*this); }
};


template<typename Hierarchy, typename Model, typename Enable = void>
class ModelView;


template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, std::enable_if_t<solver::is_hybrid_model_v<Model>>>
    : public BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>
{
    using Super        = BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>;
    using VecField     = typename Model::vecfield_type;
    using TensorFieldT = Model::ions_type::tensorfield_type;

public:
    using Model_t = Model;

    ModelView(Hierarchy& hierarchy, Model& model)
        : Super{hierarchy, model}
    {
        declareMomentumTensorAlgos();
    }

    NO_DISCARD VecField& getB() const { return this->model_.state.electromag.B; }

    NO_DISCARD VecField& getE() const { return this->model_.state.electromag.E; }

    NO_DISCARD auto& getIons() const { return this->model_.state.ions; }

    void fillPopMomTensor(auto& lvl, auto const time, auto const popidx)
    {
        using value_type = TensorFieldT::value_type;
        auto constexpr N = core::detail::tensor_field_dim_from_rank<2>();

        auto& rm   = *(this->model_.resourcesManager);
        auto& ions = this->model_.state.ions;

        for (auto patch : rm.enumerate(lvl, ions, sumTensor_))
            for (std::uint8_t c = 0; c < N; ++c)
                std::memcpy(sumTensor_[c].data(), ions[popidx].momentumTensor()[c].data(),
                            ions[popidx].momentumTensor()[c].size() * sizeof(value_type));

        MTAlgos[popidx].getOrCreateSchedule(this->hierarchy_, lvl.getLevelNumber()).fillData(time);

        for (auto patch : rm.enumerate(lvl, ions, sumTensor_))
            for (std::uint8_t c = 0; c < N; ++c)
                std::memcpy(ions[popidx].momentumTensor()[c].data(), sumTensor_[c].data(),
                            ions[popidx].momentumTensor()[c].size() * sizeof(value_type));
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(); }

protected:
    void declareMomentumTensorAlgos()
    {
        auto& rm = *(this->model_.resourcesManager);

        auto const dst_name = sumTensor_.name();

        for (auto& pop : this->model_.state.ions)
        {
            auto& MTAlgo        = MTAlgos.emplace_back();
            auto const src_name = pop.momentumTensor().name();

            auto&& [idDst, idSrc] = rm.getIDsList(dst_name, src_name);
            MTAlgo.MTalgo->registerRefine(
                idDst, idSrc, idDst, nullptr,
                std::make_shared<
                    amr::TensorFieldGhostInterpOverlapFillPattern<typename Super::GridLayoutT,
                                                                  /*rank_=*/2>>());
        }

        // can't create schedules here as the hierarchy has no levels yet
    }

    struct MTAlgo
    {
        auto& getOrCreateSchedule(auto& hierarchy, int const ilvl)
        {
            if (not MTschedules.count(ilvl))
                MTschedules.try_emplace(
                    ilvl,
                    MTalgo->createSchedule(hierarchy.getPatchLevel(ilvl), 0,
                                           std::make_shared<amr::FieldBorderSumTransactionFactory<
                                               typename Super::TensorFieldData_t>>()));
            return *MTschedules[ilvl];
        }

        std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> MTalgo
            = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> MTschedules;
    };

    std::vector<MTAlgo> MTAlgos;
    TensorFieldT sumTensor_{"PHARE_sumTensor", core::HybridQuantity::Tensor::M};
};


template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, std::enable_if_t<solver::is_mhd_model_v<Model>>>
    : public BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>
{
    using Field    = typename Model::field_type;
    using VecField = typename Model::vecfield_type;

public:
    using Model_t = Model;
    using BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>::BaseModelView;

    NO_DISCARD const Field& getRho() const { return this->model_.state.rho; }

    NO_DISCARD const VecField& getRhoV() const { return this->model_.state.rhoV; }

    NO_DISCARD const VecField& getB() const { return this->model_.state.B; }

    NO_DISCARD const Field& getEtot() const { return this->model_.state.Etot; }

    NO_DISCARD const VecField& getE() const
    {
        throw std::runtime_error("E not currently available in MHD diagnostics");
    }

    // for setBuffer function in visitHierarchy
    NO_DISCARD Field& getRho() { return this->model_.state.rho; }

    NO_DISCARD VecField& getRhoV() { return this->model_.state.rhoV; }

    NO_DISCARD VecField& getB() { return this->model_.state.B; }

    NO_DISCARD Field& getEtot() { return this->model_.state.Etot; }

    NO_DISCARD VecField& getE()
    {
        throw std::runtime_error("E not currently available in MHD diagnostics");
    }

    // diag only
    NO_DISCARD VecField& getV() { return V_diag_; }

    NO_DISCARD const VecField& getV() const { return V_diag_; }

    NO_DISCARD Field& getP() { return P_diag_; }

    NO_DISCARD const Field& getP() const { return P_diag_; }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(V_diag_, P_diag_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(V_diag_, P_diag_);
    }

protected:
    // these quantities are not always up to date in the calculations but we can compute them from
    // the conservative variables when needed their registration and allocation are handled in the
    // model
    VecField V_diag_{"diagnostics_V_", core::MHDQuantity::Vector::V};
    Field P_diag_{"diagnostics_P_", core::MHDQuantity::Scalar::P};
};


} // namespace PHARE::diagnostic



#endif // DIAGNOSTIC_MODEL_VIEW_HPP
