#ifndef PHARE_AMR_PHYSICAL_MODELS_MODEL_TEMPORARIES_HPP
#define PHARE_AMR_PHYSICAL_MODELS_MODEL_TEMPORARIES_HPP

#include "core/def.hpp"
#include "core/utilities/variants.hpp"

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
 * runtime resources. They are requested and found by physical quantity, matched on its centering
 * (GridLayout::centering), never by name, so any temporary laid out as the quantity is can be used.
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

    // FieldTmp for a Scalar, else the first of TensorTmps for the Vector or Tensor quantity type
    template<typename Qty>
    using tmp_for_ = std::conditional_t<std::is_same_v<Qty, typename GridLayout::Scalar>, FieldTmp,
                                        typename tensor_tmp_for<Qty, TensorTmps...>::type>;

public:
    using Resources = std::variant<FieldTmp, TensorTmps...>;

    NO_DISCARD auto& getRunTimeResourcesViewList() { return resources; }
    NO_DISCARD auto& getRunTimeResourcesViewList() const { return resources; }


    // Scalar for a field, Vector for a vecfield, Tensor for a tensorfield, matched on the
    // centering of qty
    template<typename Qty>
        requires(not std::is_void_v<tmp_for_<Qty>>)
    NO_DISCARD auto& getTmp(Qty const qty)
    {
        return getTmp_<tmp_for_<Qty>>(GridLayout::centering(qty));
    }

    // the next temporary of the same type and centering as previous
    template<typename Tmp>
    NO_DISCARD Tmp& getNextTmp(Tmp const& previous)
    {
        return getTmp_<Tmp>(GridLayout::centering(previous),
                            core::index_in_variants(resources, previous) + 1);
    }


    // Scalar, Vector or Tensor quantities of the model, see getTmp and the class description
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

    // takes the first temporary of the kind not yet used by this call, or creates it for qty
    template<typename Qty>
        requires(not std::is_void_v<tmp_for_<Qty>>)
    void request_(Qty const qty, std::vector<std::size_t>& used)
    {
        using Tmp            = tmp_for_<Qty>;
        auto const centering = GridLayout::centering(qty);
        std::size_t nOfType  = 0;
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

        auto& tmp = std::get<Tmp>(resources.emplace_back(Tmp{tmpName_<Tmp>(nOfType), qty}));
        static_cast<Model&>(*this).resourcesManager->registerResources(tmp);
        used.push_back(resources.size() - 1);
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
