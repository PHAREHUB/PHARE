#ifndef PHARE_CONCRETE_TAGGER_HPP
#define PHARE_CONCRETE_TAGGER_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/types.hpp"

#include "tagger.hpp"
#include "tagging_criteria.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/types/amr_types.hpp"

#include "initializer/data_provider.hpp"

#include <SAMRAI/pdat/CellData.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>


namespace PHARE::amr
{
namespace tagger_detail
{
    inline std::string nameTail(std::string const& s)
    {
        auto const p = s.rfind('_');
        return p == std::string::npos ? s : s.substr(p + 1);
    }

    // "EM_B_x" -> "Bx" ; otherwise empty
    inline std::string compactComponentName(std::string const& s)
    {
        auto const p = s.rfind('_');
        if (p == std::string::npos || p == 0)
            return {};
        auto const dir = s.substr(p + 1);
        if (dir.size() != 1 || (dir[0] != 'x' && dir[0] != 'y' && dir[0] != 'z'))
            return {};
        return nameTail(s.substr(0, p)) + dir;
    }

    // Recursively walks the model's resource tree (tuples returned by
    // getCompileTimeResourcesViewList) and pushes every scalar Field whose name matches one of
    // the requested strings. A TensorField (e.g. VecField) matching by name expands to all its
    // components. A scalar Field matches by full name, last token, or compact "Bx" form.
    template<typename FieldPtr, typename Node>
    void collectFields(std::vector<FieldPtr>& out, std::vector<std::string> const& want, Node& node)
    {
        auto const matches = [&](std::string const& full, bool acceptCompact) {
            for (auto const& w : want)
            {
                if (w == full || w == nameTail(full))
                    return true;
                if (acceptCompact && w == compactComponentName(full))
                    return true;
            }
            return false;
        };

        if constexpr (requires {
                          node.components();
                          node.name();
                      })
        {
            // TensorField (rank >= 1): match by vector name, else recurse into components.
            if (matches(node.name(), /*acceptCompact=*/false))
            {
                std::apply([&](auto&... c) { (out.push_back(&c), ...); }, node.components());
            }
            else
            {
                std::apply([&](auto&... c) { (collectFields<FieldPtr>(out, want, c), ...); },
                           node.components());
            }
        }
        else if constexpr (requires {
                               node.name();
                               node.physicalQuantity();
                           })
        {
            // Scalar Field.
            if (matches(node.name(), /*acceptCompact=*/true))
                out.push_back(&node);
        }
        else if constexpr (requires { node.getCompileTimeResourcesViewList(); })
        {
            std::apply([&](auto&... c) { (collectFields<FieldPtr>(out, want, c), ...); },
                       node.getCompileTimeResourcesViewList());
        }
        // else: unrelated node (e.g. ion population vector entries) — ignored.
    }

} // namespace tagger_detail


//-----------------------------------------------------------------------------
//  ConcreteTaggerKernel: the pure criterion + boundary-mask logic for a single
//  patch. Non-polymorphic and free of SAMRAI plumbing, so it is directly
//  unit-testable against lightweight model fixtures (which only need to expose
//  get_B()/get_B1(), the optional state tree, and the inner-boundary geometry).
//-----------------------------------------------------------------------------

template<typename Model>
class ConcreteTaggerKernel
{
    using gridlayout_type           = typename Model::gridlayout_type;
    static auto constexpr dimension = Model::dimension;

    // Default magnetic field: hybrid stores B in electromag, MHD stores it directly in the
    // state; lightweight test models expose get_B(). Also seeds the field_type alias below.
    static auto& defaultB_(Model& model)
    {
        if constexpr (requires { model.get_B(); })
            return model.get_B();
        else if constexpr (solver::is_hybrid_model_v<Model>)
            return model.state.electromag.B;
        else if constexpr (solver::is_mhd_model_v<Model>)
            return model.state.B;
        else
            static_assert(core::dependent_false_v<Model>);
    }

public:
    using vecfield_t = std::decay_t<decltype(defaultB_(std::declval<Model&>()))>;
    using field_type = typename vecfield_t::field_type;

private:
    struct QuantityTag
    {
        std::string name;
        double threshold;
    };

public:
    ConcreteTaggerKernel(initializer::PHAREDict const& dict, int maxLevelNumber = 1)
        : method_{parseTaggingMethod(dict["method"].template to<std::string>())}
        , finestLevel_{maxLevelNumber - 1}
    {
        auto const nbrQuantities = dict["nbr_quantities"].template to<int>();
        for (int i = 0; i < nbrQuantities; ++i)
        {
            auto const path        = "Q" + std::to_string(i);
            std::string const name = dict[path]["name"].template to<std::string>();
            double const threshold = dict[path]["threshold"].template to<double>();
            quantities_.push_back({name, threshold});
        }

        // method-specific params, all optional (Python's check_tagging validates the keys)
        if (dict.contains("params"))
        {
            auto const& params = dict["params"];
            if (params.contains("reltol"))
                lohnerReltol_ = params["reltol"].template to<double>();
            if (params.contains("abstol"))
                lohnerAbstol_ = params["abstol"].template to<double>();
            if (params.contains("level_scaling"))
                levelScaling_ = params["level_scaling"].template to<double>() != 0.;
        }

        // physical domain volume |Omega|, used by the wavelet level-scaled threshold;
        // injected by Python (populateDict), absent in unit-test dicts (then 1).
        if (dict.contains("domain_volume"))
            domainVolume_ = dict["domain_volume"].template to<double>();
    }

    // criterion on a single patch (writes the SAMRAI tag buffer, no ghosts).
    void tagFields(Model& model, gridlayout_type const& layout, int* tags) const
    {
        tagCells_(model, layout, tags);
    }

private:
    void tagCells_(Model& model, gridlayout_type const& layout, int* tags) const;
    std::vector<field_type const*> components_(Model& model, std::string const& qty) const;

    static std::vector<field_type const*> defaultBComponents_(Model& model);

    TaggingMethod method_;
    int finestLevel_;
    std::vector<QuantityTag> quantities_;
    double lohnerReltol_ = 0.02;
    double lohnerAbstol_ = 1e-30;
    bool levelScaling_   = true;
    double domainVolume_ = 1.0;
};


//-----------------------------------------------------------------------------
//  ConcreteTagger: the Tagger the integrator registers. Owns a kernel and adds
//  the SAMRAI patch plumbing + model.tags persistence around it.
//-----------------------------------------------------------------------------

template<typename Model>
class ConcreteTagger : public Tagger
{
    using patch_t         = typename Tagger::patch_t;
    using amr_t           = PHARE::amr::SAMRAI_Types;
    using IPhysicalModel  = PHARE::solver::IPhysicalModel<amr_t>;
    using gridlayout_type = typename Model::gridlayout_type;

public:
    ConcreteTagger(initializer::PHAREDict const& dict, int maxLevelNumber = 1)
        : Tagger{Model::model_name == "HybridModel" ? "HybridTagger" : "MHDTagger"}
        , kernel_{dict, maxLevelNumber}
    {
    }

    void tag(IPhysicalModel& model, patch_t& patch, int tag_index) override;

private:
    ConcreteTaggerKernel<Model> kernel_;
};


//-----------------------------------------------------------------------------
//                           Definitions
//-----------------------------------------------------------------------------


template<typename Model>
std::vector<typename ConcreteTaggerKernel<Model>::field_type const*>
ConcreteTaggerKernel<Model>::defaultBComponents_(Model& model)
{
    using Comp = PHARE::core::Component;
    auto& B    = defaultB_(model);
    return {&B.getComponent(Comp::X), &B.getComponent(Comp::Y), &B.getComponent(Comp::Z)};
}


template<typename Model>
std::vector<typename ConcreteTaggerKernel<Model>::field_type const*>
ConcreteTaggerKernel<Model>::components_(Model& model, std::string const& qty) const
{
    // "B" keeps the legacy default semantics (B1-preferred magnetic components). Any other
    // name is resolved generically against the model's resource tree (e.g. "rho", "V", "Bx").
    if (qty == "B")
        return defaultBComponents_(model);

    std::vector<field_type const*> out;
    if constexpr (requires { model.state; })
        tagger_detail::collectFields<field_type const*>(out, {qty}, model.state);
    else
        tagger_detail::collectFields<field_type const*>(out, {qty}, model);

    if (out.empty())
        throw std::runtime_error("tagging quantity '" + qty
                                 + "' matched no field for this model");
    return out;
}


template<typename Model>
void ConcreteTaggerKernel<Model>::tagCells_(Model& model, gridlayout_type const& layout,
                                            int* tags) const
{
    auto const nbrCells = layout.nbrCells();

    // SAMRAI tags int* buffer is FORTRAN ordering so we set false to the view; it
    // has no ghost cells (one entry per physical cell).
    bool constexpr c_ordering = false;
    auto tagsv                = core::NdArrayView<dimension, int, c_ordering>(tags, nbrCells);

    // union over quantities: a cell is tagged if ANY quantity's indicator exceeds
    // its own threshold. Zero once, then only ever set to 1.
    std::fill(tags, tags + core::product(nbrCells), 0);

    // We loop on dual-cell indexes; each component is projected onto the cell center
    // by its CellCenteredSampler (see tagging_criteria.hpp). Primal and dual share the
    // same physicalStartIndex, so a primal field is read at the same integer base.
    std::array<std::uint32_t, dimension> start;
    start[0] = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::X);
    if constexpr (dimension > 1)
        start[1] = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::Y);
    if constexpr (dimension > 2)
        start[2] = layout.physicalStartIndex(core::QtyCentering::dual, core::Direction::Z);

    // The criterion reaches +/-stencilReach(method_) and the cell-center projection adds
    // +1 upward in primal directions, so a primal-in-d component touches [i-reach, i+reach+1]
    // in d. Skip the physical cells whose stencil would leave the allocated ghost box:
    // shaveLo at the bottom (needs reach ghosts), shaveHi at the top (needs reach+1).
    int constexpr nghost = static_cast<int>(gridlayout_type::nbrGhosts());
    int const reach      = stencilReach(method_);
    auto const shaveLo   = static_cast<std::uint32_t>(std::max(0, reach - nghost));
    auto const shaveHi   = static_cast<std::uint32_t>(std::max(0, reach + 1 - nghost));

    std::array<std::uint32_t, dimension> loopLo;
    std::array<std::uint32_t, dimension> loopHi; // [loopLo, loopHi) in 0-based cell index
    for (std::size_t d = 0; d < dimension; ++d)
    {
        loopLo[d] = shaveLo;
        loopHi[d] = nbrCells[d] > shaveHi ? nbrCells[d] - shaveHi : 0;
    }

    // level-scaled threshold (wavelet only): Harten's strategy (Domingues et al. 2019,
    // Eq. 7) eps_l = eps / |O| * 2^{dim (l - L)} with |O| the PHYSICAL DOMAIN volume
    // (constant across levels) and L the finest level, so refinement is triggered more
    // eagerly on coarse levels (controls the L1 norm of the discarded details).
    double thresholdScale = 1.0;
    if (method_ == TaggingMethod::Wavelet and levelScaling_)
        thresholdScale = std::pow(2.0, static_cast<int>(dimension)
                                           * (layout.levelNumber() - finestLevel_))
                         / domainVolume_;

    // wavelet sibling pairing must follow the GLOBAL (AMR) grid: parity of cell c in
    // direction d is (AMRBox.lower[d] + c) & 1.
    auto const& amrLower = layout.AMRBox().lower;

    using sampler_t = CellCenteredSampler<dimension, field_type>;

    for (auto const& q : quantities_)
    {
        auto const comps = components_(model, q.name);

        // one cell-center sampler per component, carrying that component's own centering.
        std::vector<sampler_t> samplers;
        samplers.reserve(comps.size());
        for (auto const* c : comps)
        {
            auto const ctr = layout.centering(*c);
            std::array<bool, dimension> isPrimal;
            for (std::size_t d = 0; d < dimension; ++d)
                isPrimal[d] = (ctr[d] == core::QtyCentering::primal);
            samplers.emplace_back(*c, isPrimal);
        }
        std::vector<sampler_t const*> sampPtrs;
        sampPtrs.reserve(samplers.size());
        for (auto const& s : samplers)
            sampPtrs.push_back(&s);

        // `cell` is the 0-based physical cell index (also the tag buffer index)
        auto const indicator = [&](std::array<std::uint32_t, dimension> const& cell) {
            std::array<std::uint32_t, dimension> idx;
            for (std::size_t d = 0; d < dimension; ++d)
                idx[d] = start[d] + cell[d];

            switch (method_)
            {
                case TaggingMethod::Lohner:
                    return lohnerIndicator<dimension>(sampPtrs, idx, lohnerReltol_,
                                                      lohnerAbstol_);
                case TaggingMethod::Wavelet: {
                    std::array<std::uint32_t, dimension> parity;
                    for (std::size_t d = 0; d < dimension; ++d)
                        parity[d] = static_cast<std::uint32_t>(
                            (amrLower[d] + static_cast<int>(cell[d])) & 1);
                    return waveletIndicator<dimension>(sampPtrs, idx, parity);
                }
                default: return defaultIndicator<dimension>(sampPtrs, idx);
            }
        };

        auto const threshold = q.threshold * thresholdScale;

        if constexpr (dimension == 1)
        {
            for (std::uint32_t iCell = loopLo[0]; iCell < loopHi[0]; ++iCell)
                if (indicator({iCell}) > threshold)
                    tagsv(iCell) = 1;
        }
        else if constexpr (dimension == 2)
        {
            for (std::uint32_t ix = loopLo[0]; ix < loopHi[0]; ++ix)
                for (std::uint32_t iy = loopLo[1]; iy < loopHi[1]; ++iy)
                    if (indicator({ix, iy}) > threshold)
                        tagsv(ix, iy) = 1;
        }
        else if constexpr (dimension == 3)
        {
            for (std::uint32_t ix = loopLo[0]; ix < loopHi[0]; ++ix)
                for (std::uint32_t iy = loopLo[1]; iy < loopHi[1]; ++iy)
                    for (std::uint32_t iz = loopLo[2]; iz < loopHi[2]; ++iz)
                        if (indicator({ix, iy, iz}) > threshold)
                            tagsv(ix, iy, iz) = 1;
        }
    }
}


template<typename Model>
void ConcreteTagger<Model>::tag(PHARE::solver::IPhysicalModel<amr_t>& model, patch_t& patch,
                                int tag_index)
{
    auto& concreteModel = dynamic_cast<Model&>(model);
    auto layout         = PHARE::amr::layoutFromPatch<gridlayout_type>(patch);
    auto modelIsOnPatch = concreteModel.setOnPatch(patch);
    auto pd   = dynamic_cast<SAMRAI::pdat::CellData<int>*>(patch.getPatchData(tag_index).get());
    auto tags = pd->getPointer();

    kernel_.tagFields(concreteModel, layout, tags);


    // These tags will be saved even if they are not used in diags during this advance
    // concreteModel.tags may contain vectors for patches and levels that no longer exist
    auto key = std::to_string(patch.getPatchLevelNumber()) + "_"
               + amr::to_string(patch.getGlobalId());

    auto nCells = core::product(layout.nbrCells());

    bool item_exists_and_valid
        = concreteModel.tags.count(key) and concreteModel.tags[key]->size() == nCells;

    if (!item_exists_and_valid)
    {
        using Map_value_type = typename std::decay_t<decltype(concreteModel.tags)>::mapped_type;


        concreteModel.tags[key]
            = std::make_shared<typename Map_value_type::element_type>(layout.nbrCells());
    }

    auto nbrCells = layout.nbrCells();
    auto tagsv    = core::NdArrayView<Model::dimension, int>(concreteModel.tags[key]->data(),
                                                             layout.nbrCells());
    auto tagsvF   = core::NdArrayView<Model::dimension, int, false>(tags, layout.nbrCells());
    if constexpr (Model::dimension == 2)
    {
        for (auto iTag_x = 0u; iTag_x < nbrCells[0]; ++iTag_x)
        {
            for (auto iTag_y = 0u; iTag_y < nbrCells[1]; ++iTag_y)
            {
                tagsv(iTag_x, iTag_y) = tagsvF(iTag_x, iTag_y);
            }
        }
    }
}

} // namespace PHARE::amr

#endif // PHARE_CONCRETE_TAGGER_HPP
