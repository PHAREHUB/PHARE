#ifndef DIAGNOSTIC_MODEL_VIEW_HPP
#define DIAGNOSTIC_MODEL_VIEW_HPP

#include "core/def.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "mpi/mpi_utils.hpp"

#include "amr/amr_constants.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/messengers/field_operate_transaction.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"

#include "dict.hpp"

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
    using ResMan            = Model::resources_manager_type;
    using Field             = Model::field_type;
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
    void onLevels(Action&& action, std::size_t const minlvl = 0,
                  std::size_t const maxlvl = amr::MAX_LEVEL_IDX)
    {
        amr::onLevels(hierarchy_, std::forward<Action>(action), minlvl, maxlvl);
    }


    template<typename OnLevel, typename OrMissing>
    void onLevels(OnLevel&& onLevel, OrMissing&& orMissing, std::size_t const minlvl,
                  std::size_t const maxlvl)
    {
        amr::onLevels(hierarchy_, std::forward<OnLevel>(onLevel),
                      std::forward<OrMissing>(orMissing), minlvl, maxlvl);
    }


    template<typename Action>
    void visitHierarchy(Action&& action, int minLevel = 0, int maxLevel = 0)
    {
        amr::visitHierarchy<GridLayout>(hierarchy_, *model_.resourcesManager,
                                        std::forward<Action>(action), minLevel, maxLevel, *this,
                                        model_);
    }

    NO_DISCARD auto boundaryConditions() const { return hierarchy_.boundaryConditions(); }
    NO_DISCARD auto domainBox() const { return hierarchy_.domainBox(); }
    NO_DISCARD auto origin() const { return std::vector<double>(dimension, 0); }
    NO_DISCARD auto cellWidth() const { return hierarchy_.cellWidth(); }
    NO_DISCARD auto maxLevel() const { return hierarchy_.maxLevel(); }

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
        dict["mpi_rank"] = static_cast<std::size_t>(mpi::rank());
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

    auto operator()() const { return model_.getCompileTimeResourcesViewList(); }


    auto& field_reducer(auto& f)
        requires(core::is_field_tile_set_v<Field>)
    {
        return core::reduce_single(Model::tmpField, f);
    }
    auto& field_reducer(auto& f) { return f; }

    auto& vec_field_reducer(auto& tf_in)
        requires(core::is_field_tile_set_v<Field>)
    {
        for (std::size_t i = 0; i < 3; ++i)
            core::reduce_single(Model::tmpVec[i], tf_in[i]);
        return Model::tmpVec;
    }
    auto& vec_field_reducer(auto& f) { return f; }

    auto& tensor_field_reducer(auto& f)
        requires(core::is_field_tile_set_v<Field>)
    {
        for (std::size_t i = 0; i < 6; ++i)
            core::reduce_single(Model::tmpTensor[i], f[i]);
        return Model::tmpTensor;
    }
    auto& tensor_field_reducer(auto& f) { return f; }

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
    using Field        = Model::field_type;
    using VecField     = Model::vecfield_type;
    using TensorFieldT = Model::ions_type::tensorfield_type;

    using ScratchTensorFieldData_t
        = amr::TensorFieldData<2, typename Super::GridLayout,
                               typename Model::Super::storage_t::grid_type, core::HybridQuantity>;

public:
    using Model_t                = Model;
    using physical_quantity_type = Model::physical_quantity_type;

    ModelView(Hierarchy& hierarchy, Model& model)
        : Super{hierarchy, model}
    {
        declareMomentumTensorAlgos();
    }

    NO_DISCARD VecField& getB() const { return this->model_.state.electromag.B; }

    NO_DISCARD VecField& getE() const { return this->model_.state.electromag.E; }

    NO_DISCARD auto& getIons() const { return this->model_.state.ions; }

    auto& tmpField() { return Model::tmpField; }

    auto& tmpVecField() { return Model::tmpVec; }

    template<std::size_t rank = 2>
    auto& tmpTensorField()
    {
        static_assert(rank > 0 and rank < 3);
        if constexpr (rank == 1)
            return Model::tmpVec;
        else
            return Model::tmpTensor;
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(Model::tmpField, Model::tmpVec, Model::tmpTensor);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(Model::tmpField, Model::tmpVec, Model::tmpTensor);
    }

    void fillPopMomTensor(auto& lvl, auto const time, auto const popidx)
    {
        auto constexpr N = core::detail::tensor_field_dim_from_rank<2>();

        auto& rm   = *this->model_.resourcesManager;
        auto& ions = this->model_.state.ions;

        for (auto patch : rm.enumerate(lvl, ions, Model::tmpTensor))
            for (std::uint8_t c = 0; c < N; ++c)
                core::reduce_into(Model::tmpTensor[c], ions[popidx].momentumTensor()[c]);

        MTAlgos[popidx].getOrCreateSchedule(this->hierarchy_, lvl.getLevelNumber()).fillData(time);

        for (auto patch : rm.enumerate(lvl, ions, Model::tmpTensor))
            for (std::uint8_t c = 0; c < N; ++c)
                core::copy_fields(ions[popidx].momentumTensor()[c], Model::tmpTensor[c]);
    }


protected:
    void declareMomentumTensorAlgos()
    {
        auto& rm = *this->model_.resourcesManager;

        auto const dst_name = Model::tmpTensor.name();

        for (auto& pop : this->model_.state.ions)
        {
            auto& MTAlgo        = MTAlgos.emplace_back();
            auto const src_name = pop.momentumTensor().name();

            auto&& [idDst, idSrc] = rm.getIDsList(dst_name, src_name);
            MTAlgo.MTalgo->registerRefine(
                idDst, idSrc, idDst, nullptr,
                std::make_shared<
                    amr::TensorFieldGhostInterpOverlapFillPattern<typename Super::GridLayout,
                                                                  /*rank_=*/2>>());
        }

        // can't create schedules here as the hierarchy has no levels yet
    }

    struct MTAlgo
    {
        auto& getOrCreateSchedule(auto& hierarchy, int const ilvl)
        {
            using BorderSumOp = core::FieldBorderSumOp<typename VecField::value_type>;

            auto const level   = hierarchy.getPatchLevel(ilvl);
            auto schedule_iter = MTschedules.find(ilvl);
            auto const create_schedule
                = schedule_iter == MTschedules.end() or schedule_iter->second.level.lock() != level;

            if (create_schedule)
                schedule_iter
                    = MTschedules
                          .insert_or_assign(
                              ilvl,
                              Entry{level,
                                    MTalgo->createSchedule(
                                        level, 0,
                                        std::make_shared<amr::FieldBorderOpTransactionFactory<
                                            typename Super::TensorFieldData_t, BorderSumOp,
                                            ScratchTensorFieldData_t>>())})
                          .first;

            return *schedule_iter->second.schedule;
        }

        struct Entry
        {
            std::weak_ptr<SAMRAI::hier::PatchLevel> level; // invalidadted if the Level is destroyed
            std::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule;
        };

        std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> MTalgo
            = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();
        std::map<int, Entry> MTschedules;
    };

    std::vector<MTAlgo> MTAlgos;
};


template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, std::enable_if_t<solver::is_mhd_model_v<Model>>>
    : public BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>
{
    using Field    = Model::field_type;
    using VecField = Model::vecfield_type;

public:
    using Model_t                = Model;
    using physical_quantity_type = Model::physical_quantity_type;
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
        return std::forward_as_tuple(V_diag_, P_diag_, tmpField_, tmpVec_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(V_diag_, P_diag_, tmpField_, tmpVec_);
    }

    auto& tmpField() { return tmpField_; }

    auto& tmpVecField() { return tmpVec_; }

    template<std::size_t rank = 2>
    auto& tmpTensorField()
    {
        static_assert(rank == 1);
        return tmpVec_;
    }

protected:
    // not always current; computed from conservative vars when needed; registered/allocated in
    // model
    VecField V_diag_{"diagnostics_V_", core::MHDQuantity::Vector::V};
    Field P_diag_{"diagnostics_P_", core::MHDQuantity::Scalar::P};

    // Field/VecField may be tiled; this is reduce_single()'s destination so it must always be
    // Model's plain (never tiled) storage type, same as Model::tmpField/tmpVec.
    Model::Field_rt tmpField_{"PHARE_sumField_MHD", core::MHDQuantity::Scalar::ScalarAllPrimal};
    Model::VecField_rt tmpVec_{"PHARE_sumVec_MHD", core::MHDQuantity::Vector::VecAllPrimal};
};


} // namespace PHARE::diagnostic



#endif // DIAGNOSTIC_MODEL_VIEW_HPP
