
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_span.hpp"
#include "core/data/field/field_box_span.hpp"

#include "phare_core.hpp"
#include "phare_simulator_options.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include <chrono>
#include <cstdint>

#define PRINT(x) std::cout << __FILE__ << ":" << __LINE__ << " " << x << std::endl;

namespace PHARE::bench::core
{
static constexpr PHARE::SimOpts opts{3, 1};
using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
using GridLayout_t = TestGridLayout<typename PHARE_Types::GridLayout_t>;
using Grid_t       = PHARE_Types::Grid_t;

std::uint64_t static now()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

struct ScopeTimer
{
    ~ScopeTimer() { PRINT("Time: " << (now() - start) / div / 1e6 << " ms"); }

    double div                = 1;
    std::uint64_t const start = now();
};

} // namespace PHARE::bench::core
