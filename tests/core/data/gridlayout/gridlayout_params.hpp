#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_PARAMS_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_PARAMS_HPP

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <type_traits>


#include "core/data/grid/gridlayout.hpp"
#include "gridlayout_base_params.hpp"
#include "gridlayout_utilities.hpp"
#include "core/utilities/point/point.hpp"



using namespace PHARE::core;

/* ::std::ostream& operator<<(::std::ostream& os, GridLayoutTestParam<Layout::Yee, 1> const& param);
 */
/* void PrintTo(GridLayoutTestParam<Layout::Yee, 1> const& param, ::std::ostream& os); */


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



template<typename GridLayoutImpl>
auto createParam(std::array<floater_t<4>, GridLayoutImpl::dimension> const& dxdydz,
                 std::array<std::uint32_t, GridLayoutImpl::dimension> const& nbCellXYZ,
                 Point<floater_t<4>, GridLayoutImpl::dimension> const& origin)
{
    GridLayoutTestParam<GridLayoutImpl> param{};
    param.layout    = std::make_shared<GridLayout<GridLayoutImpl>>(dxdydz, nbCellXYZ, origin);
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
