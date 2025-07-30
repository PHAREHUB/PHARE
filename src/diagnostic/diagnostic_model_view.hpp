#ifndef DIAGNOSTIC_MODEL_VIEW_HPP
#define DIAGNOSTIC_MODEL_VIEW_HPP

#include "core/def.hpp"
#include "core/utilities/mpi_utils.hpp"

#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/messengers/field_sum_transaction.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"

#include "cppdict/include/dict.hpp"

#include <SAMRAI/xfer/RefineAlgorithm.h>

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
    using GridLayout        = Model::gridlayout_type;
    using VecField          = Model::vecfield_type;
    using TensorFieldT      = Model::ions_type::tensorfield_type;
    using GridLayoutT       = Model::gridlayout_type;
    using ResMan            = Model::resources_manager_type;
    using TensorFieldData_t = ResMan::template UserTensorField_t</*rank=*/2>::patch_data_type;
    static constexpr auto dimension = Model::dimension;


public:
    using PatchProperties
        = cppdict::Dict<float, double, std::size_t, std::vector<int>, std::vector<std::uint32_t>,
                        std::vector<double>, std::vector<std::size_t>, std::string,
                        std::vector<std::string>>;

    BaseModelView(Hierarchy& hierarchy, Model& model)
        : model_{model}
        , hierarchy_{hierarchy}
    {
        declareMomentumTensorAlgos();
    }

    NO_DISCARD std::vector<VecField*> getElectromagFields() const
    {
        return {&model_.state.electromag.B, &model_.state.electromag.E};
    }

    NO_DISCARD auto& getIons() const { return model_.state.ions; }

    void fillPopMomTensor(auto& lvl, auto const time, auto const popidx)
    {
        using value_type = TensorFieldT::value_type;
        auto constexpr N = core::detail::tensor_field_dim_from_rank<2>();

        auto& rm   = *model_.resourcesManager;
        auto& ions = model_.state.ions;

        for (auto patch : rm.enumerate(lvl, ions, sumTensor_))
            for (std::uint8_t c = 0; c < N; ++c)
                std::memcpy(sumTensor_[c].data(), ions[popidx].momentumTensor()[c].data(),
                            ions[popidx].momentumTensor()[c].size() * sizeof(value_type));

        MTAlgos[popidx].getOrCreateSchedule(hierarchy_, lvl.getLevelNumber()).fillData(time);

        for (auto patch : rm.enumerate(lvl, ions, sumTensor_))
            for (std::uint8_t c = 0; c < N; ++c)
                std::memcpy(ions[popidx].momentumTensor()[c].data(), sumTensor_[c].data(),
                            ions[popidx].momentumTensor()[c].size() * sizeof(value_type));
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

    void declareMomentumTensorAlgos()
    {
        auto& rm = *model_.resourcesManager;

        auto const dst_name = sumTensor_.name();

        for (auto& pop : model_.state.ions)
        {
            auto& MTAlgo        = MTAlgos.emplace_back();
            auto const src_name = pop.momentumTensor().name();

            auto&& [idDst, idSrc] = rm.getIDsList(dst_name, src_name);
            MTAlgo.MTalgo->registerRefine(
                idDst, idSrc, idDst, nullptr,
                std::make_shared<
                    amr::TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/2>>());
        }

        // can't create schedules here as the hierarchy has no levels yet
    }

    struct MTAlgo
    {
        auto& getOrCreateSchedule(auto& hierarchy, int const ilvl)
        {
            if (not MTschedules.count(ilvl))
                MTschedules.try_emplace(
                    ilvl, MTalgo->createSchedule(
                              hierarchy.getPatchLevel(ilvl), 0,
                              std::make_shared<
                                  amr::FieldBorderSumTransactionFactory<TensorFieldData_t>>()));
            return *MTschedules[ilvl];
        }

        std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> MTalgo
            = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> MTschedules;
    };

    std::vector<MTAlgo> MTAlgos;
    TensorFieldT sumTensor_{"PHARE_sumTensor", core::HybridQuantity::Tensor::M};
};


template<typename Hierarchy, typename Model, typename Enable = void>
class ModelView;


template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, std::enable_if_t<solver::is_hybrid_model_v<Model>>>
    : public BaseModelView<Hierarchy, Model>
{
    using VecField = typename Model::vecfield_type;

public:
    using Model_t = Model;
    using BaseModelView<Hierarchy, Model>::BaseModelView;

    NO_DISCARD std::vector<VecField*> getElectromagFields() const
    {
        return {&this->model_.state.electromag.B, &this->model_.state.electromag.E};
    }

    NO_DISCARD auto& getIons() const { return this->model_.state.ions; }
};


template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, std::enable_if_t<solver::is_mhd_model_v<Model>>>
    : public BaseModelView<Hierarchy, Model>
{
    using Field    = typename Model::field_type;
    using VecField = typename Model::vecfield_type;

public:
    using Model_t = Model;
    using BaseModelView<Hierarchy, Model>::BaseModelView;
};


} // namespace PHARE::diagnostic



#endif // DIAGNOSTIC_MODEL_VIEW_HPP
