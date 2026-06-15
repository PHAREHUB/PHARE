#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_PARAMS_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_PARAMS_HPP


#include "core/data/grid/gridlayout.hpp"
#include "core/utilities/point/point.hpp"
#include "core/models/options/hybrid_options.hpp"

#include "phare_simulator_options.hpp"

#include "gridlayout_base_params.hpp"


#include <memory>
#include <string>
#include <cassert>
#include <fstream>

using namespace PHARE::core;

/* ::std::ostream& operator<<(::std::ostream& os, GridLayoutTestParam<Layout::Yee, 1> const& param);
 */
/* void PrintTo(GridLayoutTestParam<Layout::Yee, 1> const& param, ::std::ostream& os); */


template<PHARE::SimOpts opts>
struct TestParam
{
    PHARE::HybridFieldOptions<opts> constexpr static field_options{};

    using GridLayout_t = GridLayout<PHARE::HybridOptions<field_options>{}>;
};

inline HybridQuantity::Scalar getQuantity(std::uint32_t iQuantity)
{
    switch (iQuantity)
    {
        case 0: return HybridQuantity::Scalar::Bx;
        case 1: return HybridQuantity::Scalar::By;
        case 2: return HybridQuantity::Scalar::Bz;
        case 3: return HybridQuantity::Scalar::Ex;
        case 4: return HybridQuantity::Scalar::Ey;
        case 5: return HybridQuantity::Scalar::Ez;
        case 6: return HybridQuantity::Scalar::Jx;
        case 7: return HybridQuantity::Scalar::Jy;
        case 8: return HybridQuantity::Scalar::Jz;
        case 9: return HybridQuantity::Scalar::rho;
        case 10: return HybridQuantity::Scalar::Vx;
        case 11: return HybridQuantity::Scalar::Vy;
        case 12: return HybridQuantity::Scalar::Vz;
        case 13: return HybridQuantity::Scalar::P;
        default:
            throw std::runtime_error("There is no quantity indexed by "
                                     + std::to_string(iQuantity));
    }
}



template<auto options>
auto createParam(std::array<double, options.dimension> const& dxdydz,
                 std::array<std::uint32_t, options.dimension> const& nbCellXYZ,
                 Point<double, options.dimension> const& origin)
{
    using GridLayout_t = TestParam<options>::GridLayout_t;
    GridLayoutTestParam<options> param{};
    param.layout    = std::make_shared<GridLayout_t>(dxdydz, nbCellXYZ, origin);
    param.dxdydz    = dxdydz;
    param.nbCellXYZ = nbCellXYZ;
    param.origin    = origin;
    return param;
}


template<typename Array>
bool writeToArray(std::ifstream& stream, Array& array)
{
    std::size_t dim = array.size();
    for (std::size_t i = 0; i < dim; ++i)
        if (!(stream >> array[i]))
            return false;
    return true;
}


#endif
