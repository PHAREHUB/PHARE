#ifndef PHARE_AMR_PHYSICAL_MODELS_MODEL_TEMPORARIES_HPP
#define PHARE_AMR_PHYSICAL_MODELS_MODEL_TEMPORARIES_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/variants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include <array>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <variant>
#include <type_traits>

namespace PHARE::solver
{
template<typename Tmp, typename Qty>
concept IsTensorTmpFor
    = requires { typename Tmp::tensor_t; } and std::is_same_v<typename Tmp::tensor_t, Qty>;

// the first of Tmps whose tensor quantity type is Qty (Vector: rank 1, Tensor: rank 2), else void
template<typename Qty, typename... Tmps>
struct tensor_tmp_for
{
    using type = void;
};

template<typename Qty, typename Tmp, typename... Tmps>
struct tensor_tmp_for<Qty, Tmp, Tmps...>
{
    using type = std::conditional_t<IsTensorTmpFor<Tmp, Qty>, Tmp,
                                    typename tensor_tmp_for<Qty, Tmps...>::type>;
};


/**
 * @brief ModelTemporaries holds the temporaries of a model, registered and allocated with it as
 * runtime resources. They are found by centering, never by name or physical quantity, so any
 * temporary laid out as asked for can be used.
 *
 * Users of temporaries (messengers, diagnostics) request what they need with
 * requestTemporaries() before the model is first allocated, nothing exists unless requested.
 * Within one call, repeated requests of the same kind (type and centering) are distinct
 * temporaries. Across calls they are shared: a kind only gets more temporaries when a single call
 * asks for more of it than exist, so users must not hold a temporary between calls into each
 * other.
 *
 * Model is the deriving model (CRTP), FieldTmp the field temporary type, TensorTmps the
 * vecfield/tensorfield temporary types.
 */
template<typename Model, typename GridLayout, typename FieldTmp, typename... TensorTmps>
class ModelTemporaries
{
    static_assert(not requires { typename FieldTmp::tensor_t; }, "first Tmp must be the field");

    static constexpr auto dimension = GridLayout::dimension;

public:
    using Resources = std::variant<FieldTmp, TensorTmps...>;

    NO_DISCARD auto& getRunTimeResourcesViewList() { return resources; }
    NO_DISCARD auto& getRunTimeResourcesViewList() const { return resources; }


    NO_DISCARD auto& getTmp(core::QtyCentering const centering)
    {
        return getTmp(core::ConstArray<core::QtyCentering, dimension>(centering));
    }

    NO_DISCARD auto& getTmp(std::array<core::QtyCentering, dimension> const& centering)
    {
        return getTmp_<FieldTmp>(centering);
    }

    // Vector for a vecfield, Tensor for a tensorfield, matched on the centering of qty
    template<typename Qty>
        requires(not std::is_void_v<typename tensor_tmp_for<Qty, TensorTmps...>::type>)
    NO_DISCARD auto& getTmp(Qty const qty)
    {
        using Tmp = tensor_tmp_for<Qty, TensorTmps...>::type;
        return getTmp_<Tmp>(GridLayout::centering(qty));
    }

    // the next temporary of the same type and centering as previous
    template<typename Tmp>
    NO_DISCARD Tmp& getNextTmp(Tmp const& previous)
    {
        return getTmp_<Tmp>(GridLayout::centering(previous),
                            core::index_in_variants(resources, previous) + 1);
    }


    // QtyCentering or std::array<QtyCentering, dim> for a field, Vector or Tensor quantities for a
    //  vecfield or tensorfield, see class description
    template<typename... Requests>
    void requestTemporaries(Requests const&... requests)
    {
        std::vector<std::size_t> used; // indices taken by this call
        (request_(requests, used), ...);
    }

protected:
    std::vector<Resources> resources;
    bool temporariesAllocated_ = false; // set by Model::allocate, no new temporaries after

private:
    template<typename Tmp>
    Tmp& getTmp_(auto const& centering, std::size_t const start = 0)
    {
        return core::get_from_variants_if<Tmp>(
            resources, [&](auto const& tmp) { return GridLayout::centering(tmp) == centering; },
            start);
    }

    void request_(core::QtyCentering const centering, std::vector<std::size_t>& used)
    {
        request_(core::ConstArray<core::QtyCentering, dimension>(centering), used);
    }

    void request_(std::array<core::QtyCentering, dimension> const& centering,
                  std::vector<std::size_t>& used)
    {
        ensure_<FieldTmp>(centering, used, [&]() { return fieldQuantityFor_(centering); });
    }

    template<typename Qty>
        requires(not std::is_void_v<typename tensor_tmp_for<Qty, TensorTmps...>::type>)
    void request_(Qty const qty, std::vector<std::size_t>& used)
    {
        using Tmp = tensor_tmp_for<Qty, TensorTmps...>::type;
        ensure_<Tmp>(GridLayout::centering(qty), used, [&]() { return qty; });
    }

    // takes the first temporary of the kind not yet used by this call, or creates it
    template<typename Tmp>
    void ensure_(auto const& centering, std::vector<std::size_t>& used, auto&& quantity)
    {
        std::size_t nOfType = 0;
        for (std::size_t idx = 0; idx < resources.size(); ++idx)
            if (auto const* tmp = std::get_if<Tmp>(&resources[idx]))
            {
                ++nOfType;
                if (GridLayout::centering(*tmp) == centering
                    and std::find(used.begin(), used.end(), idx) == used.end())
                {
                    used.push_back(idx);
                    return;
                }
            }

        if (temporariesAllocated_)
            throw std::runtime_error("requestTemporaries: cannot create a temporary after the "
                                     "model is allocated, request before");

        auto& tmp = std::get<Tmp>(resources.emplace_back(Tmp{tmpName_<Tmp>(nOfType), quantity()}));
        static_cast<Model&>(*this).resourcesManager->registerResources(tmp);
        used.push_back(resources.size() - 1);
    }

    // any scalar quantity with the centering, the temporary is only ever found by centering
    auto fieldQuantityFor_(std::array<core::QtyCentering, dimension> const& centering) const
    {
        using PhysicalQuantity = Model::physical_quantity_type;
        using Scalar           = PhysicalQuantity::Scalar;

        if (centering
            == core::ConstArray<core::QtyCentering, dimension>(core::QtyCentering::primal))
            return PhysicalQuantity::all_primal_field;

        for (std::size_t i = 0; i < static_cast<std::size_t>(Scalar::count); ++i)
            if (GridLayout::centering(static_cast<Scalar>(i)) == centering)
                return static_cast<Scalar>(i);

        throw std::runtime_error("requestTemporaries: no scalar quantity with this centering");
    }

    // prefixed by model, SAMRAI variable names are shared across ResourcesManagers
    template<typename Tmp>
    NO_DISCARD static std::string tmpName_(std::size_t const idx)
    {
        auto const type = [] {
            if constexpr (std::is_same_v<Tmp, FieldTmp>)
                return "Field";
            else if constexpr (Tmp::rank == 1)
                return "VecField";
            else
                return "TensorField";
        }();
        return Model::model_name + "_tmp" + type + std::to_string(idx);
    }
};

} // namespace PHARE::solver

#endif /* PHARE_AMR_PHYSICAL_MODELS_MODEL_TEMPORARIES_HPP */
