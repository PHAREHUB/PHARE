#ifndef PHARE_LOHNER_TAGGER_STRATEGY_HPP
#define PHARE_LOHNER_TAGGER_STRATEGY_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "initializer/data_provider.hpp"

#include "tagger_strategy.hpp"
#include "tagging_quantities.hpp"

#include <array>
#include <cmath>
#include <tuple>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <algorithm>

namespace PHARE::amr
{
namespace lohner
{
    template<std::size_t dim, typename Component>
    double cellCentered(Component const& component, std::array<std::uint32_t, dim> const& idx)
    {
        auto const& centering = component.centering();

        double value = 0.;
        for (unsigned corner = 0; corner < (1u << dim); ++corner)
        {
            auto node     = idx;
            double weight = 1.;
            bool inCell   = true;
            for (std::size_t d = 0; d < dim; ++d)
            {
                bool const upper = (corner >> d) & 1u;
                if (centering[d] == core::QtyCentering::primal)
                {
                    weight *= 0.5;
                    node[d] += upper;
                }
                else if (upper)
                {
                    inCell = false;
                    break;
                }
            }
            if (inCell)
                value += weight * std::apply([&](auto... i) { return component(i...); }, node);
        }
        return value;
    }


    template<std::size_t dim, typename Component>
    double indicator(std::vector<Component> const& components,
                     std::array<std::uint32_t, dim> const& idx, double reltol, double abstol)
    {
        double crit = 0.;
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto im = idx;
            auto ip = idx;
            im[d] -= 1;
            ip[d] += 1;

            double num2 = 0., d1sq = 0., fltsq = 0.;
            for (auto const& component : components)
            {
                auto const a_m = cellCentered<dim>(component, im);
                auto const a_0 = cellCentered<dim>(component, idx);
                auto const a_p = cellCentered<dim>(component, ip);
                auto const num = a_p - 2.0 * a_0 + a_m;
                auto const d1  = std::abs(a_p - a_0) + std::abs(a_0 - a_m);
                auto const flt = std::abs(a_m) + 2.0 * std::abs(a_0) + std::abs(a_p);
                num2 += num * num;
                d1sq += d1 * d1;
                fltsq += flt * flt;
            }
            crit = std::max(crit, std::sqrt(num2)
                                      / (std::sqrt(d1sq) + reltol * std::sqrt(fltsq) + abstol));
        }
        return crit;
    }
} // namespace lohner


template<typename Model>
class LohnerTaggerStrategy : public TaggerStrategy<Model>
{
    using gridlayout_type           = Model::gridlayout_type;
    static auto constexpr dimension = Model::dimension;

public:
    LohnerTaggerStrategy(initializer::PHAREDict const& dict)
        : quantities_{tagging::parseQuantities(dict)}
        , reltol_{cppdict::get_value(dict, "params/reltol", 0.02)}
        , abstol_{cppdict::get_value(dict, "params/abstol", 1e-30)}
    {
    }

    void tag(Model& model, gridlayout_type const& layout, int* tags) const override;

private:
    std::vector<tagging::QuantityThreshold> quantities_;
    double reltol_;
    double abstol_;
};


template<typename Model>
void LohnerTaggerStrategy<Model>::tag(Model& model, gridlayout_type const& layout, int* tags) const
{
    auto const quantities = tagging::resolveQuantities(model, quantities_);
    auto const nbrCells   = layout.nbrCells();

    std::array<std::uint32_t, dimension> start;
    for (std::size_t d = 0; d < dimension; ++d)
        start[d]
            = layout.physicalStartIndex(core::QtyCentering::dual, static_cast<core::Direction>(d));

    bool constexpr c_ordering = false;
    auto tagsv                = core::NdArrayView<dimension, int, c_ordering>(tags, nbrCells);

    auto const tagCell = [&](auto... cell) {
        std::array<std::uint32_t, dimension> idx{cell...};
        for (std::size_t d = 0; d < dimension; ++d)
            idx[d] += start[d];

        auto const exceeds = [&](auto const& quantity) {
            return lohner::indicator<dimension>(quantity.first, idx, reltol_, abstol_)
                   > quantity.second;
        };
        tagsv(cell...) = std::any_of(std::begin(quantities), std::end(quantities), exceeds) ? 1 : 0;
    };

    if constexpr (dimension == 1)
        for (std::uint32_t ix = 0; ix < nbrCells[0]; ++ix)
            tagCell(ix);
    if constexpr (dimension == 2)
        for (std::uint32_t ix = 0; ix < nbrCells[0]; ++ix)
            for (std::uint32_t iy = 0; iy < nbrCells[1]; ++iy)
                tagCell(ix, iy);
    if constexpr (dimension == 3)
        for (std::uint32_t ix = 0; ix < nbrCells[0]; ++ix)
            for (std::uint32_t iy = 0; iy < nbrCells[1]; ++iy)
                for (std::uint32_t iz = 0; iz < nbrCells[2]; ++iz)
                    tagCell(ix, iy, iz);
}

} // namespace PHARE::amr

#endif
