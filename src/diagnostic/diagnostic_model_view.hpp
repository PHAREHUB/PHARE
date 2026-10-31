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


template<typename Hierarchy, typename Model>
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
                                        std::forward<Action>(action), minLevel, maxLevel, model_);
    }

    NO_DISCARD auto periodicities() const
    {
        auto const& periodicities = hierarchy_.periodicities();
        return std::vector<int>(periodicities.begin(), periodicities.end());
    }
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

    // the reducers' temporaries for xs, nothing if untiled
    void requestReducers(auto const&... xs)
    {
        if constexpr (core::is_field_tile_set_v<Field>)
            model_.requestTemporaries(reducer_kind_(xs)...);
    }

    // requested in one call with the reducer's: distinct when both are primal, see primal_dst
    void requestReducedToPrimal(auto const& x, auto const primal_kind)
    {
        if constexpr (core::is_field_tile_set_v<Field>)
            model_.requestTemporaries(reducer_kind_(x), primal_kind);
        else
            model_.requestTemporaries(primal_kind);
    }

    // where to convert reduced to primal into, never reduced itself
    auto& primal_dst(auto const& reduced, auto& primal_tmp)
    {
        if (static_cast<void const*>(&reduced) == static_cast<void const*>(&primal_tmp))
            return model_.getNextTmp(primal_tmp);
        return primal_tmp;
    }

protected:
    // tiled fields are reduced into an untiled temporary of the same centering
    auto& reduce_(auto& x)
    {
        using X = std::decay_t<decltype(x)>;
        if constexpr (is_tensor_<X>)
        {
            if constexpr (core::is_field_tile_set_v<typename X::field_type>)
            {
                auto& tmp = model_.getTmp(x.physicalQuantity());
                for (std::size_t i = 0; i < x.size(); ++i)
                    core::reduce_single(tmp[i], x[i]);
                return tmp;
            }
            else
                return x;
        }
        else if constexpr (core::is_field_tile_set_v<X>)
            return core::reduce_single(model_.getTmp(GridLayout::centering(x.physicalQuantity())),
                                       x);
        else
            return x;
    }

    Model& model_;
    Hierarchy& hierarchy_;

private:
    template<typename X>
    static constexpr bool is_tensor_ = requires { typename X::tensor_t; };

    auto reducer_kind_(auto const& x) const
    {
        if constexpr (is_tensor_<std::decay_t<decltype(x)>>)
            return x.physicalQuantity();
        else
            return GridLayout::centering(x.physicalQuantity());
    }
};


template<typename Hierarchy, typename Model, typename Enable = void>
class ModelView;


template<typename Hierarchy, typename Model>
class ModelView<Hierarchy, Model, std::enable_if_t<solver::is_hybrid_model_v<Model>>>
    : public BaseModelView<Hierarchy, Model>
{
    using Super        = BaseModelView<Hierarchy, Model>;
    using Field        = Model::field_type;
    using VecField     = Model::vecfield_type;
    using TensorFieldT = Model::ions_type::tensorfield_type;

    using ScratchTensorFieldData_t
        = amr::TensorFieldData<2, typename Super::GridLayout, typename Model::storage_t::grid_type,
                               core::HybridQuantity>;

public:
    using Model_t                = Model;
    using physical_quantity_type = Model::physical_quantity_type;

    ModelView(Hierarchy& hierarchy, Model& model)
        : Super{hierarchy, model}
    {
    }

    NO_DISCARD VecField& getB() const { return this->model_.state.electromag.B; }

    NO_DISCARD VecField& getE() const { return this->model_.state.electromag.E; }

    NO_DISCARD auto& getIons() const { return this->model_.state.ions; }

    auto& tmpField() { return this->model_.getTmp(core::QtyCentering::primal); }

    auto& tmpVecField() { return this->model_.getTmp(core::HybridQuantity::Vector::V); }

    template<std::size_t rank = 2>
    auto& tmpTensorField()
    {
        static_assert(rank > 0 and rank < 3);
        if constexpr (rank == 1)
            return tmpVecField();
        else
            return this->model_.getTmp(core::HybridQuantity::Tensor::M);
    }

    // requests for the temporaries above, see TypeWriter::requestTemporaries
    void requestTmpField() { this->model_.requestTemporaries(core::QtyCentering::primal); }

    void requestTmpVecField() { this->model_.requestTemporaries(core::HybridQuantity::Vector::V); }

    // x reduced (if tiled) then converted to primal, see BaseModelView::primal_dst
    void requestTmpFieldFor(auto const& x)
    {
        this->requestReducedToPrimal(x, core::QtyCentering::primal);
    }
    void requestTmpVecFieldFor(auto const& x)
    {
        this->requestReducedToPrimal(x, core::HybridQuantity::Vector::V);
    }

    // scratch for fillPopMomTensor
    void requestMomentumTensor()
    {
        this->model_.requestTemporaries(core::HybridQuantity::Tensor::M);
        if (MTAlgos.empty())
            declareMomentumTensorAlgos();
    }

    void fillPopMomTensor(auto& lvl, auto const time, auto const popidx)
    {
        auto constexpr N = core::detail::tensor_field_dim_from_rank<2>();

        auto& rm   = *this->model_.resourcesManager;
        auto& ions = this->model_.state.ions;
        auto& tmp  = tmpTensorField();

        for (auto patch : rm.enumerate(lvl, ions, tmp))
            for (std::uint8_t c = 0; c < N; ++c)
                core::reduce_into(tmp[c], ions[popidx].momentumTensor()[c]);

        MTAlgos[popidx].getOrCreateSchedule(this->hierarchy_, lvl.getLevelNumber()).fillData(time);

        for (auto patch : rm.enumerate(lvl, ions, tmp))
            for (std::uint8_t c = 0; c < N; ++c)
                core::copy_fields(ions[popidx].momentumTensor()[c], tmp[c]);
    }


    // tiled fields are reduced into the model's untiled temporaries for writing
    auto& field_reducer(auto& f) { return this->reduce_(f); }
    auto& vec_field_reducer(auto& f) { return this->reduce_(f); }
    auto& tensor_field_reducer(auto& f) { return this->reduce_(f); }

protected:
    void declareMomentumTensorAlgos()
    {
        auto& rm = *this->model_.resourcesManager;

        auto const dst_name = tmpTensorField().name();

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
                              ilvl, Entry{level,
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
    : public BaseModelView<Hierarchy, Model>
{
    using Field    = Model::field_type;
    using VecField = Model::vecfield_type;

public:
    using Model_t                = Model;
    using physical_quantity_type = Model::physical_quantity_type;
    using BaseModelView<Hierarchy, Model>::BaseModelView;

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

    // the writers compute V and P from the conservative variables into the state's primitive
    //  ones (onto tiles if tiled, reduced for writing like the hybrid momentum tensor), these are
    //  recomputed from the conservative variables at every step
    NO_DISCARD VecField& getV() { return this->model_.state.V; }

    NO_DISCARD Field& getP() { return this->model_.state.P; }

    // requests for the reducers of getV/getP and below, see TypeWriter::requestTemporaries
    void requestV() { this->requestReducers(getV()); }

    void requestP() { this->requestReducers(getP()); }

    void requestTmpField() { this->model_.requestTemporaries(core::QtyCentering::primal); }

    void requestTmpVecField()
    {
        this->model_.requestTemporaries(core::MHDQuantity::Vector::VecAllPrimal);
    }

    // x reduced (if tiled) then converted to primal, see BaseModelView::primal_dst
    void requestTmpFieldFor(auto const& x)
    {
        this->requestReducedToPrimal(x, core::QtyCentering::primal);
    }
    void requestTmpVecFieldFor(auto const& x)
    {
        this->requestReducedToPrimal(x, core::MHDQuantity::Vector::VecAllPrimal);
    }

    auto& tmpField() { return this->model_.getTmp(core::QtyCentering::primal); }

    auto& tmpVecField() { return this->model_.getTmp(core::MHDQuantity::Vector::VecAllPrimal); }

    template<std::size_t rank = 2>
    auto& tmpTensorField()
    {
        static_assert(rank == 1);
        return tmpVecField();
    }

    // tiled fields are reduced into the model's untiled temporaries for writing
    auto& field_reducer(auto& f) { return this->reduce_(f); }
    auto& vec_field_reducer(auto& f) { return this->reduce_(f); }
};


} // namespace PHARE::diagnostic



#endif // DIAGNOSTIC_MODEL_VIEW_HPP
