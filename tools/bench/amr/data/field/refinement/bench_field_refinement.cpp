

#include "core/utilities/types.hpp"

#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/magnetic_field_init_refiner.hpp"

#include "phare_core.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include <chrono>
#include <limits>
#include <cstdint>

#define PRINT(x) std::cout << __FILE__ << ":" << __LINE__ << " " << x << std::endl;

namespace PHARE::amr::bench
{

constexpr static auto ndim = 3;
constexpr static auto NaN  = std::numeric_limits<double>::quiet_NaN();
constexpr static PHARE::SimOpts opts{ndim, 1};

using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
using GridLayout_t = TestGridLayout<typename PHARE_Types::GridLayout_t>;
using Grid_t       = PHARE_Types::Grid_t;
using Refiner      = MagneticFieldInitRefiner<ndim>;


std::uint64_t static now()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

struct ScopeTimer
{
    ~ScopeTimer() { PRINT("Time: " << (now() - start) / 1e6 << " ms"); }

    double div                = 1;
    std::uint64_t const start = now();
};

int test()
{
    std::size_t constexpr static FOR = 100;

    GridLayout_t coarse_layout{.1, 300ul};
    GridLayout_t fine_layout = GridLayout_t::make(
        Box<int, ndim>{core::ConstArray<int, ndim>(10), core::ConstArray<int, ndim>(590)}, .1 / 2);

    SAMRAI::hier::Box const coarse_box = samrai_box_from(coarse_layout.AMRBox());
    SAMRAI::hier::Box const fine_box   = samrai_box_from(fine_layout.AMRBox());
    SAMRAI::hier::IntVector const ratio{SAMRAI::tbox::Dimension{ndim}, 2};

    Grid_t dst{"rho", fine_layout, PHARE::core::HybridQuantity::Scalar::rho, 1};
    Grid_t const src{"rho", coarse_layout, PHARE::core::HybridQuantity::Scalar::rho, 1};
    auto const& qty = dst.physicalQuantity();

    {
        ScopeTimer timer{FOR};
        for (std::size_t L = 0; L < FOR; ++L)
        {
            Refiner refiner{fine_layout.centering(qty), fine_box, coarse_box, ratio};
            refine_field(dst, src, fine_box, refiner);
        }
    }

    auto const mid = dst.size() / 2;
    return dst.data()[mid - 1] != dst.data()[mid + 1];
}

} // namespace PHARE::amr::bench

int main()
{
    std::ios::sync_with_stdio(false);
    auto const ret = PHARE::amr::bench::test();
    PRINT(ret);
    return ret;
}
