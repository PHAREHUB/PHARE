#ifndef PHARE_CORE_NUMERICS_RIEMANN_SOLVERS_UCT_COEFFICIENTS_HPP
#define PHARE_CORE_NUMERICS_RIEMANN_SOLVERS_UCT_COEFFICIENTS_HPP

namespace PHARE::core
{
// Which of its upwind coefficients (aL, aR, dL, dR) a Riemann solver stores for the constrained
// transport. dL is always stored. When the advective ones are 1/2 on both sides (Rusanov), aL
// and aR are not stored and read as 1/2; when the diffusive ones are equal (Rusanov, HLL), dR
// is not stored and read from dL.
struct UCTCoefficients
{
    bool halfAdvective  = false;
    bool equalDiffusive = false;
};
} // namespace PHARE::core

#endif
