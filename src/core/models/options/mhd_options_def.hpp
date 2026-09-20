#ifndef PHARE_OPTIONS_MHD_OPTIONS_DEF_HPP
#define PHARE_OPTIONS_MHD_OPTIONS_DEF_HPP

#include <cstdint>

// if mhd is off, axes stay at MHDOff and no selector specializations exist for it (canary)
namespace PHARE::MHDOpts
{

enum class TimeIntegratorType : uint8_t { MHDOff, Euler, TVDRK2, TVDRK3, SSPRK4_5, count };
enum class ReconstructionType : uint8_t { MHDOff, Constant, Linear, WENO3, WENOZ, MP5, count };
enum class SlopeLimiterType : uint8_t { MHDOff, None, VanLeer, MinMod, count };
enum class RiemannSolverType : uint8_t { MHDOff, Rusanov, HLL, HLLD, count };

} // namespace PHARE::MHDOpts

#endif // PHARE_OPTIONS_MHD_OPTIONS_DEF_HPP
