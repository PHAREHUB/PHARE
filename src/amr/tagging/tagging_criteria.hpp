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
//  Components are co-indexed on the same cell index regardless of their true
//  centering (B is primal in its own direction, rho is cell-/node-centered):
//  this is the legacy behavior and is exact for cell-centered quantities; for
//  primal components it carries a harmless half-cell offset in the vector norm.
//  Tagging is only a refinement hint, so we deliberately avoid any projection.
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
