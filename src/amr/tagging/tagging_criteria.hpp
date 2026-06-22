#ifndef PHARE_TAGGING_CRITERIA_HPP
#define PHARE_TAGGING_CRITERIA_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>
#include <stdexcept>

namespace PHARE::amr
{
//-----------------------------------------------------------------------------
//  Tagging criteria (pure, unit-testable)
//
//  Both criteria use a CENTERED stencil (i-1, i, i+1): its reach is +/-1, so the
//  evaluation is valid at every physical cell for any quantity as long as
//  nbrGhosts >= 1 (always true) -- no per-quantity ghost clamp is ever needed,
//  and it reaches one fewer ghost layer than a forward (i,i+1,i+2) stencil,
//  which matters for cell-centered (ddd) fields like the MHD mass density.
//
//  The criteria are evaluated per cell and produce a cell-centered indicator.
//  Components must be sampled AT the cell center: a component is generally
//  staggered (B is primal in its own direction, hybrid rho is all-primal, MHD
//  rho is all-dual), so reading it at the dual cell index would offset it half a
//  cell in every primal direction and break the symmetry of the centered stencil.
//  CellCenteredSampler (below) projects each component onto the cell center on the
//  fly (order-2 averaging of the bracketing primal nodes, dual directions read
//  directly) with no temporary buffer. The criteria below are centering-agnostic:
//  they call elem(i...) on whatever sampler they are given.
//
//  The criteria generalize over an arbitrary list of field components
//  (B -> 3 components, rho -> 1) and combine over directions by max.
//-----------------------------------------------------------------------------

enum class TaggingMethod { Default, Lohner };

inline TaggingMethod parseTaggingMethod(std::string const& s)
{
    if (s == "default")
        return TaggingMethod::Default;
    if (s == "lohner")
        return TaggingMethod::Lohner;
    throw std::runtime_error("unknown tagging method '" + s + "' (expected 'default' or 'lohner')");
}


// idx shifted by n along direction dir (cast through int to allow negative n)
template<std::size_t dim>
std::array<std::uint32_t, dim> shiftIndex(std::array<std::uint32_t, dim> idx, std::size_t dir, int n)
{
    idx[dir] = static_cast<std::uint32_t>(static_cast<int>(idx[dir]) + n);
    return idx;
}

// f(idx...) as double, unpacking the index array over Field::operator()
template<std::size_t dim, typename Field>
double at(Field const& f, std::array<std::uint32_t, dim> const& idx)
{
    return std::apply([&](auto... i) { return static_cast<double>(f(i...)); }, idx);
}


// Samples a (possibly staggered) field AT the cell center, projecting on the fly:
// in each PRIMAL direction it averages the two bracketing nodes (order-2,
// PrimalToDual: at dual cell i -> 0.5*(F[i] + F[i+1])); DUAL directions are read
// directly. No temporary buffer. An all-dual mask is an exact passthrough.
//
// The +1 reach in primal directions stacks with the criteria's +/-1 centered
// stencil: a primal-in-d component touches [i-1, i+2] in d, so the top physical
// cell needs nbrGhosts >= 2 (the caller guards this).
template<std::size_t dim, typename Field>
class CellCenteredSampler
{
public:
    CellCenteredSampler(Field const& f, std::array<bool, dim> const& isPrimal)
        : f_{&f}
        , isPrimal_{isPrimal}
    {
    }

    template<typename... Indices>
    double operator()(Indices... is) const
    {
        static_assert(sizeof...(Indices) == dim, "index count must match dimension");
        std::array<std::uint32_t, dim> const idx{static_cast<std::uint32_t>(is)...};

        // tensor product over the bracketing nodes: dual dirs contribute {offset 0,
        // coef 1}; primal dirs contribute {offset 0 and +1, coef 0.5 each}.
        double acc = 0.0;
        for (unsigned mask = 0; mask < (1u << dim); ++mask)
        {
            double coef = 1.0;
            auto j      = idx;
            bool valid  = true;
            for (std::size_t d = 0; d < dim; ++d)
            {
                bool const bit = (mask >> d) & 1u;
                if (isPrimal_[d])
                {
                    coef *= 0.5;
                    if (bit)
                        j[d] += 1;
                }
                else if (bit)
                {
                    valid = false; // dual direction has no +1 corner
                    break;
                }
            }
            if (valid)
                acc += coef * at<dim>(*f_, j);
        }
        return acc;
    }

private:
    Field const* f_;
    std::array<bool, dim> isPrimal_;
};


// "default" criterion: centered analog of the legacy per-component smoothness
// ratio. Per component, per direction:
//   crit_d(F) = |F(i+1)-F(i-1)| / (1 + 1/2(|F(i+1)-F(i)| + |F(i)-F(i-1)|))
// On a linear ramp of slope g this equals 2g/(1+g), reproducing the legacy
// forward operator |F(i+2)-F(i)|/(1+|F(i+1)-F(i)|) exactly. Final indicator =
// max over all components and all directions.
template<std::size_t dim, typename Field>
double defaultIndicator(std::vector<Field const*> const& comps,
                        std::array<std::uint32_t, dim> const& idx)
{
    double crit = 0.0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        auto const im = shiftIndex(idx, d, -1);
        auto const ip = shiftIndex(idx, d, +1);
        for (auto const* fp : comps)
        {
            auto const& f   = *fp;
            auto const fm   = at<dim>(f, im);
            auto const f0   = at<dim>(f, idx);
            auto const fp1  = at<dim>(f, ip);
            auto const num  = std::abs(fp1 - fm);
            auto const den  = 1.0 + 0.5 * (std::abs(fp1 - f0) + std::abs(f0 - fm));
            crit            = std::max(crit, num / den);
        }
    }
    return crit;
}


// "lohner" criterion: Loehner scale-invariant second-difference estimator, the
// centered form (port of tools/tagging_study lohner, combine=max). Per
// direction d, over components: N = ||a_p - 2 a_0 + a_m||,
//   D = || |a_p-a_0| + |a_0-a_m| || + eps * || |a_m| + 2|a_0| + |a_p| || + eps_abs,
//   crit_d = N / D ; final indicator = max over directions.
template<std::size_t dim, typename Field>
double lohnerIndicator(std::vector<Field const*> const& comps,
                       std::array<std::uint32_t, dim> const& idx, double eps = 0.02,
                       double eps_abs = 1e-30)
{
    double crit = 0.0;
    for (std::size_t d = 0; d < dim; ++d)
    {
        auto const im = shiftIndex(idx, d, -1);
        auto const ip = shiftIndex(idx, d, +1);
        double num2 = 0.0, d1sq = 0.0, fltsq = 0.0;
        for (auto const* fp : comps)
        {
            auto const& f   = *fp;
            auto const a_m  = at<dim>(f, im);
            auto const a_0  = at<dim>(f, idx);
            auto const a_p  = at<dim>(f, ip);
            auto const numc = a_p - 2.0 * a_0 + a_m;
            auto const d1c  = std::abs(a_p - a_0) + std::abs(a_0 - a_m);
            auto const fltc = std::abs(a_m) + 2.0 * std::abs(a_0) + std::abs(a_p);
            num2 += numc * numc;
            d1sq += d1c * d1c;
            fltsq += fltc * fltc;
        }
        auto const val = std::sqrt(num2) / (std::sqrt(d1sq) + eps * std::sqrt(fltsq) + eps_abs);
        crit           = std::max(crit, val);
    }
    return crit;
}

} // namespace PHARE::amr

#endif // PHARE_TAGGING_CRITERIA_HPP
